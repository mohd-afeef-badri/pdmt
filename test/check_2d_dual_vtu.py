#!/usr/bin/env python3
"""Topology checks for PDMT's two-dimensional polygonal VTU output."""

import argparse
from collections import Counter
from pathlib import Path
import xml.etree.ElementTree as ET


VTK_LINE = 3
VTK_POLYGON = 7


def data_arrays(parent):
    return {
        array.get("Name"): (array.text or "").split()
        for array in parent.findall("./DataArray")
    }


def read_mesh(path):
    piece = ET.parse(path).getroot().find(".//Piece")
    if piece is None:
        raise ValueError("VTU has no UnstructuredGrid Piece")

    point_values = list(
        map(float, (piece.find("./Points/DataArray").text or "").split())
    )
    points = [
        tuple(point_values[index : index + 3])
        for index in range(0, len(point_values), 3)
    ]

    arrays = data_arrays(piece.find("./Cells"))
    connectivity = list(map(int, arrays["connectivity"]))
    offsets = list(map(int, arrays["offsets"]))
    types = list(map(int, arrays["types"]))

    cells = []
    begin = 0
    for end, cell_type in zip(offsets, types):
        cells.append((cell_type, connectivity[begin:end]))
        begin = end
    return points, cells


def edge(a, b):
    return (a, b) if a < b else (b, a)


def check_mesh(path, expected_polygons=None, expected_boundary=None):
    points, cells = read_mesh(path)
    polygons = [nodes for cell_type, nodes in cells if cell_type == VTK_POLYGON]
    lines = [nodes for cell_type, nodes in cells if cell_type == VTK_LINE]

    if expected_polygons is not None and len(polygons) != expected_polygons:
        raise AssertionError(
            f"{path}: expected {expected_polygons} polygons, found {len(polygons)}"
        )

    polygon_edges = Counter()
    total_polygon_vertices = 0
    for cell_index, polygon in enumerate(polygons):
        if len(polygon) < 3 or len(set(polygon)) != len(polygon):
            raise AssertionError(f"{path}: polygon {cell_index} is degenerate")

        twice_area = 0.0
        for index, node in enumerate(polygon):
            next_node = polygon[(index + 1) % len(polygon)]
            polygon_edges[edge(node, next_node)] += 1
            x0, y0, _ = points[node]
            x1, y1, _ = points[next_node]
            twice_area += x0 * y1 - x1 * y0
        if abs(twice_area) < 1.0e-12:
            raise AssertionError(f"{path}: polygon {cell_index} has zero area")
        total_polygon_vertices += len(polygon)

    nonmanifold = [segment for segment, count in polygon_edges.items() if count > 2]
    if nonmanifold:
        raise AssertionError(
            f"{path}: {len(nonmanifold)} polygon segments occur more than twice"
        )

    polygon_boundary = {
        segment for segment, count in polygon_edges.items() if count == 1
    }
    explicit_boundary = {
        edge(line[0], line[1])
        for line in lines
        if len(line) == 2
    }
    if polygon_boundary != explicit_boundary:
        missing = polygon_boundary - explicit_boundary
        extra = explicit_boundary - polygon_boundary
        raise AssertionError(
            f"{path}: polygon boundary differs from explicit boundary "
            f"({len(missing)} unexpected, {len(extra)} missing segments)"
        )

    if expected_boundary is not None and len(polygon_boundary) != expected_boundary:
        raise AssertionError(
            f"{path}: expected {expected_boundary} boundary segments, "
            f"found {len(polygon_boundary)}"
        )

    return len(polygons), total_polygon_vertices, len(polygon_boundary)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mesh", type=Path)
    parser.add_argument("--expect-polygons", type=int)
    parser.add_argument("--expect-boundary", type=int)
    parser.add_argument("--more-connectivity-than", type=Path)
    args = parser.parse_args()

    polygon_count, connectivity_count, boundary_count = check_mesh(
        args.mesh, args.expect_polygons, args.expect_boundary
    )

    if args.more_connectivity_than is not None:
        _, reference_connectivity, _ = check_mesh(args.more_connectivity_than)
        if connectivity_count <= reference_connectivity:
            raise AssertionError(
                f"{args.mesh}: expected more polygon connectivity than "
                f"{args.more_connectivity_than} ({connectivity_count} <= "
                f"{reference_connectivity})"
            )

    print(
        f"{args.mesh}: {polygon_count} valid polygons, "
        f"{connectivity_count} polygon vertices, "
        f"{boundary_count} boundary segments"
    )


if __name__ == "__main__":
    main()
