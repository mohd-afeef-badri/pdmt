#!/usr/bin/env python3
"""Regression checks for the labeled 2D viscous-layer fixture."""

from collections import Counter
import math
from pathlib import Path
import sys
import xml.etree.ElementTree as ET

from check_2d_dual_vtu import VTK_LINE, VTK_POLYGON, check_mesh, read_mesh


def close(first, second):
    return math.isclose(first, second, rel_tol=0.0, abs_tol=1.0e-6)


def main():
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {Path(sys.argv[0]).name} MESH.vtu")
    path = Path(sys.argv[1])

    polygon_count, _, boundary_count = check_mesh(
        path, expected_polygons=25, expected_boundary=28
    )
    points, cells = read_mesh(path)
    if len(points) != 42:
        raise AssertionError(f"{path}: expected 42 points, found {len(points)}")

    polygons = [nodes for cell_type, nodes in cells if cell_type == VTK_POLYGON]
    lines = [nodes for cell_type, nodes in cells if cell_type == VTK_LINE]
    if any(len(polygon) != 4 for polygon in polygons[-20:]):
        raise AssertionError(f"{path}: the 20 layer cells must be quadrilaterals")
    total_area = 0.0
    for polygon in polygons:
        twice_area = sum(
            points[first][0] * points[second][1]
            - points[second][0] * points[first][1]
            for first, second in zip(polygon, polygon[1:] + polygon[:1])
        )
        total_area += 0.5 * abs(twice_area)
    if not close(total_area, 1.0):
        raise AssertionError(
            f"{path}: deformed cells cover area {total_area}, expected 1"
        )

    piece = ET.parse(path).getroot().find(".//Piece")
    label_array = piece.find("./CellData/DataArray[@Name='label']")
    labels = list(map(int, map(float, (label_array.text or "").split())))
    if labels[:polygon_count] != [20] * polygon_count:
        raise AssertionError(f"{path}: viscous cells did not retain fluid label 20")
    boundary_labels = Counter(labels[polygon_count : polygon_count + len(lines)])
    if boundary_labels != Counter({11: 2, 12: 2, 13: 24}):
        raise AssertionError(
            f"{path}: exterior physical labels changed: {boundary_labels}"
        )

    expected_rows = {
        1.0 - 0.06 * layer: [0.0, 0.5, 1.0]
        for layer in range(1, 11)
    }
    for y, expected_x in expected_rows.items():
        actual_x = sorted(x for x, point_y, _ in points[12:] if close(point_y, y))
        if len(actual_x) != len(expected_x) or any(
            not close(actual, expected)
            for actual, expected in zip(actual_x, expected_x)
        ):
            raise AssertionError(
                f"{path}: expected layer row y={y} at x={expected_x}, "
                f"found {actual_x}"
            )

    print(
        f"{path}: {polygon_count} polygons, {boundary_count} labeled boundary "
        "segments, ten full-width layers with total thickness 0.6"
    )


if __name__ == "__main__":
    main()
