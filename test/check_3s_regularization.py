#!/usr/bin/env python3
"""Check that 3S dual-area regularization improves an open surface fixture."""

import argparse
from collections import Counter
import math
from pathlib import Path
import xml.etree.ElementTree as ET


VTK_POLYGON = 7


def arrays(parent):
    return {
        item.get("Name"): (item.text or "").split()
        for item in parent.findall("./DataArray")
    }


def edge(first, second):
    return (first, second) if first < second else (second, first)


def cross(first, second):
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def read_surface(path):
    piece = ET.parse(path).getroot().find(".//Piece")
    if piece is None:
        raise AssertionError(f"{path}: missing UnstructuredGrid Piece")

    point_values = list(
        map(float, (piece.find("./Points/DataArray").text or "").split())
    )
    points = [
        tuple(point_values[index : index + 3])
        for index in range(0, len(point_values), 3)
    ]
    cell_arrays = arrays(piece.find("./Cells"))
    connectivity = list(map(int, cell_arrays["connectivity"]))
    offsets = list(map(int, cell_arrays["offsets"]))
    types = list(map(int, cell_arrays["types"]))

    polygons = []
    begin = 0
    for end, cell_type in zip(offsets, types):
        if cell_type == VTK_POLYGON:
            polygons.append(connectivity[begin:end])
        begin = end
    return points, polygons


def area_statistics(path):
    points, polygons = read_surface(path)
    edge_counts = Counter(
        edge(polygon[index], polygon[(index + 1) % len(polygon)])
        for polygon in polygons
        for index in range(len(polygon))
    )

    areas = []
    boundary = []
    for polygon in polygons:
        if len(polygon) < 3 or len(set(polygon)) != len(polygon):
            raise AssertionError(f"{path}: degenerate polygon")
        area_vector = [0.0, 0.0, 0.0]
        is_boundary = False
        for index, node in enumerate(polygon):
            next_node = polygon[(index + 1) % len(polygon)]
            segment = edge(node, next_node)
            contribution = cross(points[node], points[next_node])
            for coordinate in range(3):
                area_vector[coordinate] += contribution[coordinate]
            is_boundary |= edge_counts[segment] == 1
        area = 0.5 * math.sqrt(sum(value * value for value in area_vector))
        if area <= 1.0e-14:
            raise AssertionError(f"{path}: zero-area polygon")
        areas.append(area)
        boundary.append(is_boundary)

    mean = sum(areas) / len(areas)
    cv = math.sqrt(
        sum((area - mean) ** 2 for area in areas) / len(areas)
    ) / mean
    boundary_areas = [
        area for area, is_boundary in zip(areas, boundary) if is_boundary
    ]
    interior_areas = [
        area for area, is_boundary in zip(areas, boundary) if not is_boundary
    ]
    if not boundary_areas or not interior_areas:
        raise AssertionError(f"{path}: fixture must have boundary and interior cells")
    ratio = (
        sum(boundary_areas) / len(boundary_areas)
    ) / (
        sum(interior_areas) / len(interior_areas)
    )
    return points, polygons, cv, ratio


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("unsmoothed", type=Path)
    parser.add_argument("regularized", type=Path)
    args = parser.parse_args()

    _, original_polygons, original_cv, original_ratio = area_statistics(
        args.unsmoothed
    )
    regularized_points, regularized_polygons, regularized_cv, regularized_ratio = (
        area_statistics(args.regularized)
    )
    if [len(cell) for cell in regularized_polygons] != [
        len(cell) for cell in original_polygons
    ]:
        raise AssertionError("regularization changed the polygon topology")
    if regularized_cv >= original_cv - 1.0e-6:
        raise AssertionError(
            "regularization did not reduce surface-cell area variation "
            f"({original_cv:.6g} -> {regularized_cv:.6g})"
        )
    if abs(regularized_ratio - 1.0) >= abs(original_ratio - 1.0):
        raise AssertionError(
            "regularization did not bring mean boundary area closer to "
            f"the interior area ({original_ratio:.6g} -> "
            f"{regularized_ratio:.6g})"
        )

    corners = {
        (0.0, 0.0, 0.0),
        (1.0, 0.0, 0.0),
        (1.0, 1.0, 0.0),
        (0.0, 1.0, 0.0),
    }
    for corner in corners:
        if not any(
            all(abs(point[index] - corner[index]) < 1.0e-12
                for index in range(3))
            for point in regularized_points
        ):
            raise AssertionError(
                f"regularization did not preserve boundary corner {corner}"
            )

    print(
        f"{args.regularized}: area CV {original_cv:.6f} -> "
        f"{regularized_cv:.6f}; boundary/interior ratio "
        f"{original_ratio:.6f} -> {regularized_ratio:.6f}"
    )


if __name__ == "__main__":
    main()
