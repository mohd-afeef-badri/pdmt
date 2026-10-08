#!/usr/bin/env python3
"""Check PDMT VTU face planarity, simplicity, and optional volume conservation."""

import argparse
import collections
import itertools
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def subtract(a, b):
    return tuple(a[index] - b[index] for index in range(3))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def dot(a, b):
    return sum(a[index] * b[index] for index in range(3))


def norm(vector):
    return math.sqrt(dot(vector, vector))


def read_explicit_faces(path):
    root = ET.parse(path).getroot()
    piece = root.find(".//Piece")
    if piece is None:
        raise AssertionError(f"{path}: missing UnstructuredGrid Piece")

    point_array = piece.find("./Points/DataArray")
    if point_array is None:
        raise AssertionError(f"{path}: missing point coordinates")
    values = list(map(float, (point_array.text or "").split()))
    points = [tuple(values[index:index + 3])
              for index in range(0, len(values), 3)]

    fields = {
        item.get("Name"): (item.text or "").split()
        for item in root.findall(".//FieldData/DataArray")
    }
    connectivity = list(map(int, fields["pdmt_face_connectivity"]))
    offsets = list(map(int, fields["pdmt_face_offsets"]))
    faces = []
    begin = 0
    for end in offsets:
        faces.append(connectivity[begin:end])
        begin = end
    return points, faces


def read_polyhedron_cells(path):
    root = ET.parse(path).getroot()
    arrays = {
        item.get("Name"): (item.text or "").split()
        for item in root.findall(".//Cells/DataArray")
    }
    stream = list(map(int, arrays["faces"]))
    offsets = list(map(int, arrays["faceoffsets"]))
    cells = []
    begin = 0
    for end in offsets:
        cursor = begin
        face_count = stream[cursor]
        cursor += 1
        cell = []
        for _ in range(face_count):
            node_count = stream[cursor]
            cursor += 1
            cell.append(stream[cursor:cursor + node_count])
            cursor += node_count
        if cursor != end:
            raise AssertionError("invalid VTK polyhedron face stream")
        cells.append(cell)
        begin = end
    return cells


def polyhedron_volume(points, faces):
    signed_volume = 0.0
    for face in faces:
        origin = points[face[0]]
        for index in range(1, len(face) - 1):
            signed_volume += dot(
                origin,
                cross(points[face[index]], points[face[index + 1]]),
            ) / 6.0
    return abs(signed_volume)


def normalized_warp(points, face):
    coordinates = [points[index] for index in face]
    diameter = max(
        norm(subtract(a, b))
        for a, b in itertools.combinations(coordinates, 2)
    )
    if diameter == 0.0:
        raise AssertionError("face has zero diameter")

    # Use the largest-area triangle so the reference plane remains stable for
    # polygons containing collinear circumcentric subdivision points.
    best_origin = None
    best_normal = None
    best_area2 = 0.0
    for a, b, c in itertools.combinations(coordinates, 3):
        normal = cross(subtract(b, a), subtract(c, a))
        area2 = norm(normal)
        if area2 > best_area2:
            best_origin = a
            best_normal = normal
            best_area2 = area2
    if best_area2 <= 1.0e-14 * diameter * diameter:
        raise AssertionError("face is degenerate or collinear")

    maximum_distance = max(
        abs(dot(best_normal, subtract(point, best_origin))) / best_area2
        for point in coordinates
    )
    return maximum_distance / diameter


def has_self_intersection(points, face, tolerance=1.0e-12):
    if len(face) <= 3:
        return False
    coordinates = [points[index] for index in face]
    best_normal = max(
        (cross(subtract(b, a), subtract(c, a))
         for a, b, c in itertools.combinations(coordinates, 3)),
        key=norm,
    )
    dropped = max(range(3), key=lambda component: abs(best_normal[component]))
    retained = [component for component in range(3) if component != dropped]
    projected = [(point[retained[0]], point[retained[1]])
                 for point in coordinates]
    scale = max(
        math.dist(a, b) for a, b in itertools.combinations(projected, 2)
    )
    length_tolerance = tolerance * scale
    area_tolerance = tolerance * scale * scale

    def orientation(a, b, c):
        return ((b[0] - a[0]) * (c[1] - a[1]) -
                (b[1] - a[1]) * (c[0] - a[0]))

    def on_segment(point, a, b):
        return (abs(orientation(a, b, point)) <= area_tolerance and
                min(a[0], b[0]) - length_tolerance <= point[0] <=
                max(a[0], b[0]) + length_tolerance and
                min(a[1], b[1]) - length_tolerance <= point[1] <=
                max(a[1], b[1]) + length_tolerance)

    def sign(value):
        if value > area_tolerance:
            return 1
        if value < -area_tolerance:
            return -1
        return 0

    edge_count = len(projected)
    for first in range(edge_count):
        first_next = (first + 1) % edge_count
        a, b = projected[first], projected[first_next]
        if math.dist(a, b) <= length_tolerance:
            return True
        for second in range(first + 1, edge_count):
            second_next = (second + 1) % edge_count
            if second == first_next or second_next == first:
                continue
            c, d = projected[second], projected[second_next]
            abc, abd = orientation(a, b, c), orientation(a, b, d)
            cda, cdb = orientation(c, d, a), orientation(c, d, b)
            if ((sign(abc) * sign(abd) < 0 and
                 sign(cda) * sign(cdb) < 0) or
                    (sign(abc) == 0 and on_segment(c, a, b)) or
                    (sign(abd) == 0 and on_segment(d, a, b)) or
                    (sign(cda) == 0 and on_segment(a, c, d)) or
                    (sign(cdb) == 0 and on_segment(b, c, d))):
                return True
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--tolerance", type=float, default=1.0e-11)
    parser.add_argument("--check-simple", action="store_true")
    parser.add_argument("--check-closed", action="store_true")
    parser.add_argument("--expected-volume", type=float)
    parser.add_argument("--volume-tolerance", type=float, default=1.0e-10)
    args = parser.parse_args()

    points, faces = read_explicit_faces(args.mesh)
    warps = [normalized_warp(points, face) for face in faces]
    maximum = max(warps, default=0.0)
    if maximum > args.tolerance:
        face = max(range(len(warps)), key=warps.__getitem__)
        raise AssertionError(
            f"face {face} is not planar: normalized warp {maximum:.6g} "
            f"> {args.tolerance:.6g}"
        )
    if args.check_simple:
        butterflies = [face for face in range(len(faces))
                       if has_self_intersection(points, faces[face])]
        if butterflies:
            raise AssertionError(
                f"self-intersecting polygonal faces: {butterflies[:20]}"
            )
    if args.check_closed or args.expected_volume is not None:
        cells = read_polyhedron_cells(args.mesh)
    if args.check_closed:
        for cell_id, cell in enumerate(cells):
            edge_counts = collections.Counter(
                tuple(sorted((face[index], face[(index + 1) % len(face)])))
                for face in cell
                for index in range(len(face))
            )
            invalid = [edge for edge, count in edge_counts.items()
                       if count != 2]
            if invalid:
                raise AssertionError(
                    f"polyhedron {cell_id} has an open or non-manifold shell; "
                    f"invalid edges: {invalid[:20]}"
                )
    if args.expected_volume is not None:
        volume = sum(polyhedron_volume(points, cell) for cell in cells)
        relative_error = abs(volume - args.expected_volume) / abs(
            args.expected_volume
        )
        if relative_error > args.volume_tolerance:
            raise AssertionError(
                f"dual volume {volume:.16g} does not conserve expected "
                f"volume {args.expected_volume:.16g}; relative error "
                f"{relative_error:.3e} > {args.volume_tolerance:.3e}"
            )
    print(f"{len(faces)} faces are planar; maximum normalized warp {maximum:.3e}")


if __name__ == "__main__":
    main()
