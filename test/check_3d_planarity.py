#!/usr/bin/env python3
"""Verify that every explicit polygonal face in a PDMT VTU is planar."""

import argparse
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--tolerance", type=float, default=1.0e-11)
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
    print(f"{len(faces)} faces are planar; maximum normalized warp {maximum:.3e}")


if __name__ == "__main__":
    main()
