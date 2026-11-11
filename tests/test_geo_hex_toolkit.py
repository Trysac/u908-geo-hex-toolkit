"""Tests for geo_hex_toolkit exercising public APIs, acceptance items, and oracles."""

from __future__ import annotations

import random

import h3
import pytest
import shapely
import shapely.geometry

from geo_hex_toolkit import (
    GeoHexError,
    cell_children,
    cell_parent,
    cell_resolution,
    cell_to_boundary,
    cell_to_latlng,
    cell_to_string,
    grid_disk,
    grid_distance,
    haversine_m,
    latlng_to_cell,
    polygon_to_cells,
    project,
)


def test_acceptance_latlng_to_cell_and_resolution() -> None:
    """latlng_to_cell(0.0, 0.0, 0) == 0x0800000020000000 and cell_to_string(0x0800000020000000) == '0800000020000000'; cell_resolution of a cell from latlng_to_cell(40.4168, -3.7038, 9) == 9."""
    c0 = latlng_to_cell(0.0, 0.0, 0)
    assert c0 == 0x0800000020000000
    assert cell_to_string(0x0800000020000000) == "0800000020000000"

    c9 = latlng_to_cell(40.4168, -3.7038, 9)
    assert cell_resolution(c9) == 9


def test_acceptance_cell_children_and_parent() -> None:
    """cell_children(latlng_to_cell(0.0, 0.0, 0)) returns 7 resolution-1 cells with coordinates (0, 0), (1, 0), ..."""
    c0 = latlng_to_cell(0.0, 0.0, 0)
    children_c0 = cell_children(c0)
    assert len(children_c0) == 7

    # Verify coordinates of children_c0
    expected_c0_coords = [(0, 0), (1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)]
    for child, (eq, er) in zip(children_c0, expected_c0_coords, strict=True):
        res = cell_resolution(child)
        assert res == 1
        q = ((child >> 30) & ((1 << 30) - 1)) - (1 << 29)
        r = (child & ((1 << 30) - 1)) - (1 << 29)
        assert (q, r) == (eq, er)
        assert cell_parent(child, 0) == c0

    # Cell (1, 0) at res 0: q=1, r=0
    # Center in fine coords is (2q - r, q + 3r) = (2, 1)
    # Neighbours: (2+1, 1)=(3, 1), (2+1, 1-1)=(3, 0), (2, 0), (1, 1), (1, 2), (2, 2)
    c_1_0 = (0 << 60) | ((1 + (1 << 29)) << 30) | (0 + (1 << 29))
    children_1_0 = cell_children(c_1_0)
    expected_coords = [(2, 1), (3, 1), (3, 0), (2, 0), (1, 1), (1, 2), (2, 2)]
    assert len(children_1_0) == 7
    for child, (eq, er) in zip(children_1_0, expected_coords, strict=True):
        q = ((child >> 30) & ((1 << 30) - 1)) - (1 << 29)
        r = (child & ((1 << 30) - 1)) - (1 << 29)
        assert (q, r) == (eq, er)
        assert cell_parent(child, 0) == c_1_0


def test_acceptance_grid_disk_and_grid_distance() -> None:
    """grid_disk(c, 1) has 7 cells and grid_disk(c, 2) has 19 cells for any cell c, and grid_distance(c, d) == 2 for every d in grid_disk(c, 2) not in grid_disk(c, 1)."""
    c = latlng_to_cell(48.8566, 2.3522, 5)
    d1 = grid_disk(c, 1)
    d2 = grid_disk(c, 2)
    assert len(d1) == 7
    assert len(d2) == 19
    assert d1[0] == c
    assert d2[0] == c

    set_d1 = set(d1)
    outer_ring = [cell for cell in d2 if cell not in set_d1]
    assert len(outer_ring) == 12
    for d in outer_ring:
        assert grid_distance(c, d) == 2


def test_acceptance_haversine() -> None:
    """haversine_m((0.0, 0.0), (0.0, 1.0)) == pytest.approx(111195.05197522943, rel=1e-12) and haversine_m((51.5, -0.12), (48.85, 2.35)) == pytest.approx(343128.26457181945, rel=1e-12)."""
    d1 = haversine_m((0.0, 0.0), (0.0, 1.0))
    assert d1 == pytest.approx(111195.05197522943, rel=1e-12)

    d2 = haversine_m((51.5, -0.12), (48.85, 2.35))
    assert d2 == pytest.approx(343128.26457181945, rel=1e-12)


def test_acceptance_oracles_dev() -> None:
    """Oracles (dev): for 500 seeded random points, shapely covers point; grid_disk(c, 1) shares edge; haversine equals h3."""
    rng = random.Random(42)  # noqa: S311
    points: list[tuple[float, float]] = []
    for _ in range(500):
        lat = rng.uniform(-75.0, 75.0)
        lng = rng.uniform(-170.0, 170.0)
        points.append((lat, lng))

    for lat, lng in points:
        cell = latlng_to_cell(lat, lng, 7)
        boundary = cell_to_boundary(cell)
        proj_poly = shapely.geometry.Polygon([project(v_lat, v_lng) for v_lat, v_lng in boundary])
        pt = shapely.geometry.Point(project(lat, lng))
        assert proj_poly.covers(pt)

    # Test shared edge for a cell
    test_cell = latlng_to_cell(12.34, 56.78, 6)
    poly_c = shapely.geometry.Polygon([project(la, lo) for la, lo in cell_to_boundary(test_cell)])
    # Edge length at res 6
    e6 = 1_000_000.0 / (7.0 ** (6 / 2.0))
    for neighbour in grid_disk(test_cell, 1)[1:]:
        poly_n = shapely.geometry.Polygon(
            [project(la, lo) for la, lo in cell_to_boundary(neighbour)]
        )
        inter = poly_c.intersection(poly_n)
        assert inter.length == pytest.approx(e6, rel=1e-9)

    # Test haversine_m vs h3.great_circle_distance
    for i in range(20):
        p1 = points[i]
        p2 = points[i + 20]
        dist_calc = haversine_m(p1, p2)
        dist_h3 = h3.great_circle_distance(p1, p2, unit="m")
        assert dist_calc == pytest.approx(dist_h3, rel=1e-9)


def test_acceptance_polygon_to_cells() -> None:
    """polygon_to_cells(ring, 6) returns exactly cells whose projected centres lie inside ring; raises too_many_cells."""
    ring = [
        (10.0, 10.0),
        (10.1, 10.0),
        (10.1, 10.1),
        (10.0, 10.1),
    ]
    cells = polygon_to_cells(ring, 6)
    assert len(cells) > 0
    shapely_poly = shapely.geometry.Polygon([project(lat, lng) for lat, lng in ring])
    for cell in cells:
        clat, clng = cell_to_latlng(cell)
        pt = shapely.geometry.Point(project(clat, clng))
        assert shapely_poly.contains(pt)

    # Exceeding max_cells raises GeoHexError with code 'too_many_cells'
    with pytest.raises(GeoHexError) as exc_info:
        polygon_to_cells(ring, 6, max_cells=1)
    assert exc_info.value.code == "too_many_cells"


def test_acceptance_errors() -> None:
    """latlng_to_cell(91.0, 0.0, 3) raises 'invalid_latlng'; res 16 raises 'invalid_resolution'; cell 0 raises 'invalid_cell'; 10_001 vertices raises 'input_too_large'; ring > 180 deg raises 'antimeridian_not_supported'."""
    with pytest.raises(GeoHexError) as exc_info:
        latlng_to_cell(91.0, 0.0, 3)
    assert exc_info.value.code == "invalid_latlng"

    with pytest.raises(GeoHexError) as exc_info:
        latlng_to_cell(0.0, 0.0, 16)
    assert exc_info.value.code == "invalid_resolution"

    with pytest.raises(GeoHexError) as exc_info:
        cell_to_latlng(0)
    assert exc_info.value.code == "invalid_cell"

    huge_ring = [(0.0, 0.0)] * 10_001
    with pytest.raises(GeoHexError) as exc_info:
        polygon_to_cells(huge_ring, 3)
    assert exc_info.value.code == "input_too_large"

    wide_ring = [(0.0, -100.0), (0.0, 90.0), (1.0, 0.0)]
    with pytest.raises(GeoHexError) as exc_info:
        polygon_to_cells(wide_ring, 3)
    assert exc_info.value.code == "antimeridian_not_supported"


def test_edge_cases_and_error_paths() -> None:
    """Test poles, invalid k, distance resolution mismatch, invalid parent resolution, negative axial coords, short ring."""
    # Poles
    cell_north = latlng_to_cell(90.0, 0.0, 2)
    cell_south = latlng_to_cell(-90.0, 0.0, 2)
    assert cell_resolution(cell_north) == 2
    assert cell_resolution(cell_south) == 2

    # Invalid k
    c = latlng_to_cell(0.0, 0.0, 1)
    with pytest.raises(GeoHexError) as exc:
        grid_disk(c, -1)
    assert exc.value.code == "invalid_k"

    with pytest.raises(GeoHexError) as exc:
        grid_disk(c, 101)
    assert exc.value.code == "invalid_k"

    # Distance between different resolutions
    c_res2 = latlng_to_cell(0.0, 0.0, 2)
    with pytest.raises(GeoHexError) as exc:
        grid_distance(c, c_res2)
    assert exc.value.code == "invalid_resolution"

    # Parent resolution > cell resolution
    with pytest.raises(GeoHexError) as exc:
        cell_parent(c, 2)
    assert exc.value.code == "invalid_resolution"

    # Children of resolution 15
    c15 = latlng_to_cell(0.0, 0.0, 15)
    with pytest.raises(GeoHexError) as exc:
        cell_children(c15)
    assert exc.value.code == "invalid_resolution"

    # Short ring returns empty list
    assert polygon_to_cells([(0.0, 0.0), (1.0, 1.0)], 2) == []

    # Haversine invalid latlng
    with pytest.raises(GeoHexError) as exc:
        haversine_m((95.0, 0.0), (0.0, 0.0))
    assert exc.value.code == "invalid_latlng"

    with pytest.raises(GeoHexError) as exc:
        haversine_m((0.0, 0.0), (0.0, 200.0))
    assert exc.value.code == "invalid_latlng"

    # k=0 disk returns just center
    assert grid_disk(c, 0) == [c]
