#!/usr/bin/env python3
"""Check that dual-volume regularization improves a 3D dual fixture."""

import argparse
import math
from pathlib import Path
import xml.etree.ElementTree as ET


def arrays(parent):
    return {
        item.get("Name"): (item.text or "").split()
        for item in parent.findall("./DataArray")
    }


def read_polyhedron_volumes(path):
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
    faces = list(map(int, cell_arrays["faces"]))
    face_offsets = list(map(int, cell_arrays["faceoffsets"]))

    volumes = []
    begin = 0
    for end in face_offsets:
        encoded = faces[begin:end]
        begin = end
        face_count = encoded[0]
        cursor = 1
        signed_volume = 0.0
        for _ in range(face_count):
            vertex_count = encoded[cursor]
            ids = encoded[cursor + 1 : cursor + 1 + vertex_count]
            cursor += vertex_count + 1
            a = points[ids[0]]
            for index in range(1, vertex_count - 1):
                b = points[ids[index]]
                c = points[ids[index + 1]]
                cross_bc = (
                    b[1] * c[2] - b[2] * c[1],
                    b[2] * c[0] - b[0] * c[2],
                    b[0] * c[1] - b[1] * c[0],
                )
                signed_volume += sum(a[i] * cross_bc[i] for i in range(3)) / 6.0
        volumes.append(abs(signed_volume))
    return points, volumes


def coefficient_of_variation(values):
    mean = sum(values) / len(values)
    variance = sum((value - mean) ** 2 for value in values) / len(values)
    return math.sqrt(variance) / mean


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("unsmoothed", type=Path)
    parser.add_argument("regularized", type=Path)
    args = parser.parse_args()

    _, unsmoothed_volumes = read_polyhedron_volumes(args.unsmoothed)
    regularized_points, regularized_volumes = read_polyhedron_volumes(
        args.regularized
    )
    if len(unsmoothed_volumes) != len(regularized_volumes):
        raise AssertionError("regularization changed the number of polyhedra")
    if any(volume <= 1.0e-14 for volume in regularized_volumes):
        raise AssertionError("regularization produced a zero-volume polyhedron")

    unsmoothed_cv = coefficient_of_variation(unsmoothed_volumes)
    regularized_cv = coefficient_of_variation(regularized_volumes)
    if regularized_cv >= unsmoothed_cv - 1.0e-6:
        raise AssertionError(
            "regularization did not improve cell-volume variation "
            f"({regularized_cv:.6g} >= {unsmoothed_cv:.6g})"
        )

    # The fixture has eight boundary seed cells followed by one interior cell.
    unsmoothed_ratio = (
        sum(unsmoothed_volumes[:8]) / 8.0 / unsmoothed_volumes[8]
    )
    regularized_ratio = (
        sum(regularized_volumes[:8]) / 8.0 / regularized_volumes[8]
    )
    if abs(regularized_ratio - 1.0) >= abs(unsmoothed_ratio - 1.0):
        raise AssertionError(
            "regularization did not bring mean boundary volume closer to "
            f"the interior volume ({unsmoothed_ratio:.6g} -> "
            f"{regularized_ratio:.6g})"
        )

    corner_count = sum(
        all(abs(coordinate) < 1.0e-12 or abs(coordinate - 1.0) < 1.0e-12
            for coordinate in point)
        for point in regularized_points
    )
    if corner_count != 8:
        raise AssertionError(
            f"regularization did not preserve the eight boundary corners "
            f"(found {corner_count})"
        )

    print(
        "3D regularization reduced cell-volume CV from "
        f"{unsmoothed_cv:.6f} to {regularized_cv:.6f}; "
        "mean boundary/interior ratio changed from "
        f"{unsmoothed_ratio:.6f} to {regularized_ratio:.6f}"
    )


if __name__ == "__main__":
    main()
