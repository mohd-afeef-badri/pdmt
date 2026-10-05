#!/usr/bin/env python3
"""Check the local perpendicular-bisector geometry of a 3S dual."""

import argparse
import itertools
import math
from pathlib import Path
import xml.etree.ElementTree as ET


VTK_POLYGON = 7


def subtract(first, second):
    return tuple(first[index] - second[index] for index in range(3))


def add(first, second):
    return tuple(first[index] + second[index] for index in range(3))


def scale(vector, factor):
    return tuple(value * factor for value in vector)


def dot(first, second):
    return sum(first[index] * second[index] for index in range(3))


def cross(first, second):
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def norm(vector):
    return math.sqrt(dot(vector, vector))


def edge(first, second):
    return (first, second) if first < second else (second, first)


def triangle_circumcenter(a, b, c):
    u = subtract(b, a)
    v = subtract(c, a)
    normal = cross(u, v)
    denominator = 2.0 * dot(normal, normal)
    if denominator == 0.0:
        raise AssertionError("fixture contains a degenerate triangle")
    offset = scale(
        add(
            scale(cross(v, normal), dot(u, u)),
            scale(cross(normal, u), dot(v, v)),
        ),
        1.0 / denominator,
    )
    return add(a, offset)


def read_gmsh_surface(path):
    lines = [line.strip() for line in path.read_text().splitlines()]
    node_start = lines.index("$Nodes")
    node_count = int(lines[node_start + 1])
    nodes = {}
    for line in lines[node_start + 2:node_start + 2 + node_count]:
        values = line.split()
        nodes[int(values[0])] = tuple(map(float, values[1:4]))

    element_start = lines.index("$Elements")
    element_count = int(lines[element_start + 1])
    triangles = []
    for line in lines[element_start + 2:element_start + 2 + element_count]:
        values = list(map(int, line.split()))
        if values[1] != 2:
            continue
        tag_count = values[2]
        triangles.append(tuple(values[3 + tag_count:6 + tag_count]))
    return nodes, triangles


def read_vtu(path):
    piece = ET.parse(path).getroot().find(".//Piece")
    if piece is None:
        raise AssertionError(f"{path}: missing UnstructuredGrid Piece")
    values = list(
        map(float, (piece.find("./Points/DataArray").text or "").split())
    )
    points = [
        tuple(values[index:index + 3])
        for index in range(0, len(values), 3)
    ]
    arrays = {
        item.get("Name"): (item.text or "").split()
        for item in piece.findall("./Cells/DataArray")
    }
    connectivity = list(map(int, arrays["connectivity"]))
    offsets = list(map(int, arrays["offsets"]))
    types = list(map(int, arrays["types"]))
    polygons = []
    begin = 0
    for end, cell_type in zip(offsets, types):
        if cell_type == VTK_POLYGON:
            polygons.append(connectivity[begin:end])
        begin = end
    return points, polygons


def closest_point(points, target, tolerance):
    distances = [norm(subtract(point, target)) for point in points]
    result = min(range(len(points)), key=distances.__getitem__)
    if distances[result] > tolerance:
        raise AssertionError(
            f"missing expected circumcentric point {target}; "
            f"nearest distance is {distances[result]:.6g}"
        )
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("primal", type=Path)
    parser.add_argument("dual", type=Path)
    parser.add_argument("--tolerance", type=float, default=2.0e-5)
    args = parser.parse_args()

    primal_points, triangles = read_gmsh_surface(args.primal)
    points, polygons = read_vtu(args.dual)
    if len(polygons) != len(primal_points):
        raise AssertionError("expected one closed-fan polygon per primal vertex")
    if any(len(polygon) != 6 or len(set(polygon)) != 6 for polygon in polygons):
        raise AssertionError("regular tetrahedron dual must have four hexagons")

    polygon_segments = {
        edge(polygon[index], polygon[(index + 1) % len(polygon)])
        for polygon in polygons
        for index in range(len(polygon))
    }
    maximum_orthogonality_error = 0.0
    checked_segments = 0
    for triangle in triangles:
        coordinates = [primal_points[node] for node in triangle]
        centre = triangle_circumcenter(*coordinates)
        centre_node = closest_point(points, centre, args.tolerance)
        for first, second in itertools.combinations(range(3), 2):
            midpoint = scale(add(coordinates[first], coordinates[second]), 0.5)
            midpoint_node = closest_point(points, midpoint, args.tolerance)
            if edge(centre_node, midpoint_node) not in polygon_segments:
                raise AssertionError(
                    "circumcentre-to-edge-midpoint segment is absent from dual"
                )
            primal_vector = subtract(coordinates[second], coordinates[first])
            dual_vector = subtract(points[centre_node], points[midpoint_node])
            error = abs(dot(primal_vector, dual_vector)) / (
                norm(primal_vector) * norm(dual_vector)
            )
            maximum_orthogonality_error = max(
                maximum_orthogonality_error, error
            )
            checked_segments += 1

    if maximum_orthogonality_error > args.tolerance:
        raise AssertionError(
            "surface primal/dual orthogonality error is "
            f"{maximum_orthogonality_error:.6g}"
        )
    print(
        f"{len(polygons)} valid circumcentric surface polygons; "
        f"{checked_segments} local dual segments checked; maximum "
        f"orthogonality error {maximum_orthogonality_error:.3e}"
    )


if __name__ == "__main__":
    main()
