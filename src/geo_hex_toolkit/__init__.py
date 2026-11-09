"""Pure Python aperture-7 hierarchical hexagon grid on an equal-area projection."""

from __future__ import annotations

import logging
import math
from collections.abc import Sequence
from typing import Literal

__all__ = [
    "GeoHexError",
    "GeoHexErrorCode",
    "R",
    "cell_children",
    "cell_parent",
    "cell_resolution",
    "cell_to_boundary",
    "cell_to_latlng",
    "cell_to_string",
    "grid_disk",
    "grid_distance",
    "haversine_m",
    "latlng_to_cell",
    "polygon_to_cells",
    "project",
]

logger = logging.getLogger(__name__)
logger.addHandler(logging.NullHandler())

# Authalic sphere radius (same radius value as H3 uses)
R: float = 6371007.180918475
E0: float = 1_000_000.0
SQRT3: float = math.sqrt(3.0)

# Bit layout constants
OFFSET_29: int = 1 << 29  # 2**29 = 536870912
COORD_MIN: int = -OFFSET_29
COORD_MAX: int = OFFSET_29 - 1
MASK_30: int = (1 << 30) - 1
MASK_4: int = 0xF
MASK_64: int = (1 << 64) - 1

MAX_RESOLUTION: int = 15
MAX_RING_VERTICES: int = 10_000
MAX_DISK_K: int = 100

GeoHexErrorCode = Literal[
    "invalid_latlng",
    "invalid_resolution",
    "invalid_cell",
    "invalid_k",
    "input_too_large",
    "too_many_cells",
    "antimeridian_not_supported",
]


class GeoHexError(Exception):
    """Base error class for geo_hex_toolkit operations."""

    def __init__(self, code: GeoHexErrorCode, message: str = "") -> None:
        self.code = code
        super().__init__(message or code)


# Axial neighbour offsets in specified counter-clockwise ring order:
# (1, 0), (1, -1), (0, -1), (-1, 0), (-1, 1), (0, 1)
NEIGHBOUR_OFFSETS: tuple[tuple[int, int], ...] = (
    (1, 0),
    (1, -1),
    (0, -1),
    (-1, 0),
    (-1, 1),
    (0, 1),
)


VERTEX_OFFSETS: tuple[tuple[int, int], ...] = (
    (1, 1),
    (-1, 2),
    (-2, 1),
    (-1, -1),
    (1, -2),
    (2, -1),
)


def _encode_cell(res: int, q: int, r: int) -> int:
    if not (0 <= res <= MAX_RESOLUTION):
        raise GeoHexError("invalid_resolution", f"Resolution must be 0..{MAX_RESOLUTION}")
    if not (COORD_MIN <= q <= COORD_MAX and COORD_MIN <= r <= COORD_MAX):
        raise GeoHexError("invalid_cell", "Axial coordinates out of range [-2^29, 2^29)")
    val = (res << 60) | ((q + OFFSET_29) << 30) | (r + OFFSET_29)
    return val & MASK_64


def _decode_cell(cell: int) -> tuple[int, int, int]:
    if not isinstance(cell, int) or cell <= 0 or cell > MASK_64:
        raise GeoHexError("invalid_cell", "Cell identifier must be an unsigned 64-bit integer")
    res = (cell >> 60) & MASK_4
    q_enc = (cell >> 30) & MASK_30
    r_enc = cell & MASK_30
    if res > MAX_RESOLUTION:
        raise GeoHexError("invalid_resolution", f"Invalid cell resolution {res}")
    q = q_enc - OFFSET_29
    r = r_enc - OFFSET_29
    if not (COORD_MIN <= q <= COORD_MAX and COORD_MIN <= r <= COORD_MAX):
        raise GeoHexError("invalid_cell", "Invalid decoded cell coordinates")
    return res, q, r


def cell_resolution(cell: int) -> int:
    """Return the resolution (0..15) of a cell id."""
    res, _, _ = _decode_cell(cell)
    return res


def cell_to_string(cell: int) -> str:
    """Return the 16-character hexadecimal string representation of cell."""
    _decode_cell(cell)
    return f"{cell:016x}"


def project(lat: float, lng: float) -> tuple[float, float]:
    """Project (lat, lng) to Lambert cylindrical equal-area coordinates in metres."""
    if not (-90.0 <= lat <= 90.0 and -180.0 <= lng <= 180.0):
        raise GeoHexError("invalid_latlng", f"Latitude {lat} or Longitude {lng} out of range")
    x = R * math.radians(lng)
    y = R * math.sin(math.radians(lat))
    return x, y


def _unproject(x: float, y: float) -> tuple[float, float]:
    """Convert projected metres back to (lat, lng) in degrees."""
    # Numerical clamping for pole limits
    clamped_sin = max(-1.0, min(1.0, y / R))
    lat = math.degrees(math.asin(clamped_sin))
    lng = math.degrees(x / R)
    # Normalise longitude to [-180, 180)
    lng = (lng + 180.0) % 360.0 - 180.0
    return lat, lng


def _cube_round(frac_q: float, frac_r: float) -> tuple[int, int]:
    """Deterministic cube rounding to nearest axial (q, r) coordinates."""
    frac_s = -frac_q - frac_r
    q = round(frac_q)
    r = round(frac_r)
    s = round(frac_s)

    q_diff = abs(q - frac_q)
    r_diff = abs(r - frac_r)
    s_diff = abs(s - frac_s)

    # Tie rule: drop the component with largest rounding error
    if q_diff > r_diff and q_diff > s_diff:
        q = -r - s
    elif r_diff > s_diff:
        r = -q - s
    return q, r


def _resolution_basis(res: int) -> tuple[float, float, float, float]:
    """Return the 2x2 column vectors of basis (b1_x, b2_x, b1_y, b2_y)."""
    # Matrix M = [[2, -1], [1, 3]]
    # fine basis is coarse basis multiplied by M^{-1}:
    # [b1_fine, b2_fine] = [b1_coarse, b2_coarse] * M^{-1}
    # where M^{-1} = (1/7) * [[3, 1], [-1, 2]]
    # At resolution res, fine basis = B0 * (M^{-1})^res.
    # Let Minv = [[3, 1], [-1, 2]]. Minv / 7.
    # Compute Minv^res:
    m00, m01 = 1, 0
    m10, m11 = 0, 1
    for _ in range(res):
        n00 = m00 * 3 + m01 * (-1)
        n01 = m00 * 1 + m01 * 2
        n10 = m10 * 3 + m11 * (-1)
        n11 = m10 * 1 + m11 * 2
        m00, m01, m10, m11 = n00, n01, n10, n11

    det = 7**res
    inv00, inv01 = m00 / det, m01 / det
    inv10, inv11 = m10 / det, m11 / det

    # Base basis for res 0
    # b1 = (sqrt(3) * e0, 0)
    # b2 = (sqrt(3) * e0 / 2, 1.5 * e0)
    b0_00 = SQRT3 * E0
    b0_01 = SQRT3 * E0 / 2.0
    b0_10 = 0.0
    b0_11 = 1.5 * E0

    # fine basis = B0 * M^{-res}
    b_00 = b0_00 * inv00 + b0_01 * inv10
    b_01 = b0_00 * inv01 + b0_01 * inv11
    b_10 = b0_10 * inv00 + b0_11 * inv10
    b_11 = b0_10 * inv01 + b0_11 * inv11
    return b_00, b_01, b_10, b_11


def _xy_to_qr(x: float, y: float, res: int) -> tuple[int, int]:
    # Fine basis B_res is B_0 * M^{-res}
    # x = B_res * [q, r]^T = B_0 * M^{-res} * [q, r]^T
    # So B_0^{-1} * x = M^{-res} * [q, r]^T
    # Let (q0, r0) be unrounded coordinates in res-0 lattice:
    # [q0, r0]^T = B_0^{-1} * x.
    # Since coarse = fine * M^{-1}? Wait,
    # fine = M^res * coarse_unrounded!
    b_00, b_01, b_10, b_11 = _resolution_basis(res)
    det = b_00 * b_11 - b_01 * b_10
    frac_q = (b_11 * x - b_01 * y) / det
    frac_r = (-b_10 * x + b_00 * y) / det
    # Deterministic cube rounding tie rule: drop the largest rounding error
    return _cube_round(frac_q, frac_r)


def _qr_to_xy(q: int, r: int, res: int) -> tuple[float, float]:
    b_00, b_01, b_10, b_11 = _resolution_basis(res)
    x = b_00 * q + b_01 * r
    y = b_10 * q + b_11 * r
    return x, y


def latlng_to_cell(lat: float, lng: float, res: int) -> int:
    """Index a (lat, lng) location to a hexagon cell at the given resolution (0..15)."""
    if not (0 <= res <= MAX_RESOLUTION):
        raise GeoHexError("invalid_resolution", f"Resolution must be 0..{MAX_RESOLUTION}")
    x, y = project(lat, lng)
    q, r = _xy_to_qr(x, y, res)
    return _encode_cell(res, q, r)


def cell_to_latlng(cell: int) -> tuple[float, float]:
    """Return the (lat, lng) center of the given cell."""
    res, q, r = _decode_cell(cell)
    x, y = _qr_to_xy(q, r, res)
    return _unproject(x, y)


def cell_to_boundary(cell: int) -> tuple[tuple[float, float], ...]:
    """Return 6 (lat, lng) boundary vertices in counter-clockwise order."""
    res, q, r = _decode_cell(cell)
    b_00, b_01, b_10, b_11 = _resolution_basis(res)

    vertices: list[tuple[float, float]] = []
    for dq, dr in VERTEX_OFFSETS:
        mq = 3 * q + dq
        mr = 3 * r + dr
        vx = (mq * b_00 + mr * b_01) / 3.0
        vy = (mq * b_10 + mr * b_11) / 3.0
        # Convert back using exact formula from spec:
        # lat = degrees(asin(y / R)) and lng normalised to [-180, 180)
        clamped_sin = max(-1.0, min(1.0, vy / R))
        vlat = math.degrees(math.asin(clamped_sin))
        vlng = math.degrees(vx / R)
        vlng = (vlng + 180.0) % 360.0 - 180.0
        vertices.append((vlat, vlng))
    return tuple(vertices)


def cell_parent(cell: int, res: int) -> int:
    """Return the parent cell at resolution res <= current resolution."""
    curr_res, q, r = _decode_cell(cell)
    if not (0 <= res <= curr_res):
        raise GeoHexError("invalid_resolution", f"Parent resolution {res} invalid for {curr_res}")
    curr_q, curr_r = q, r
    for _ in range(curr_res - res):
        # Parent of fine (a, b) is cube rounding of ((3a + b) / 7, (-a + 2b) / 7)
        frac_q = (3 * curr_q + curr_r) / 7.0
        frac_r = (-curr_q + 2 * curr_r) / 7.0
        curr_q, curr_r = _cube_round(frac_q, frac_r)
    return _encode_cell(res, curr_q, curr_r)


def cell_children(cell: int) -> list[int]:
    """Return the 7 children of cell at resolution res + 1."""
    res, q, r = _decode_cell(cell)
    if res >= MAX_RESOLUTION:
        raise GeoHexError("invalid_resolution", f"Cannot take children of resolution {res}")
    # Center in fine axial coordinates is (2q - r, q + 3r)
    center_q = 2 * q - r
    center_r = q + 3 * r
    target_res = res + 1

    # 7 children: center and 6 axial neighbours in counter-clockwise order
    children = [_encode_cell(target_res, center_q, center_r)]
    for dq, dr in NEIGHBOUR_OFFSETS:
        children.append(_encode_cell(target_res, center_q + dq, center_r + dr))
    return children


def grid_distance(a: int, b: int) -> int:
    """Return the grid hex distance (|dq| + |dr| + |dq + dr|) // 2 between two cells."""
    res_a, qa, ra = _decode_cell(a)
    res_b, qb, rb = _decode_cell(b)
    if res_a != res_b:
        raise GeoHexError(
            "invalid_resolution", "Cannot calculate distance at different resolutions"
        )
    dq = qa - qb
    dr = ra - rb
    return (abs(dq) + abs(dr) + abs(dq + dr)) // 2


def grid_disk(cell: int, k: int) -> list[int]:
    """Return all cells within grid distance k, center first, then rings outward."""
    if not (0 <= k <= MAX_DISK_K):
        raise GeoHexError("invalid_k", f"k must be between 0 and {MAX_DISK_K}")
    res, q, r = _decode_cell(cell)
    results = [_encode_cell(res, q, r)]
    if k == 0:
        return results

    # Generate ring by ring outward for radius 1..k
    for radius in range(1, k + 1):
        # Start at offset radius in direction 4: (-radius, radius)
        curr_q = q - radius
        curr_r = r + radius
        # Traverse 6 directions for 'radius' steps each
        for dq, dr in NEIGHBOUR_OFFSETS:
            for _ in range(radius):
                results.append(_encode_cell(res, curr_q, curr_r))
                curr_q += dq
                curr_r += dr
    return results


def haversine_m(a: tuple[float, float], b: tuple[float, float]) -> float:
    """Return great-circle distance between two (lat, lng) points in metres."""
    lat1, lng1 = a
    lat2, lng2 = b
    if not (-90.0 <= lat1 <= 90.0 and -180.0 <= lng1 <= 180.0):
        raise GeoHexError("invalid_latlng", f"Point a {a} has invalid coordinates")
    if not (-90.0 <= lat2 <= 90.0 and -180.0 <= lng2 <= 180.0):
        raise GeoHexError("invalid_latlng", f"Point b {b} has invalid coordinates")

    phi1 = math.radians(lat1)
    phi2 = math.radians(lat2)
    delta_phi = math.radians(lat2 - lat1)
    delta_lambda = math.radians(lng2 - lng1)

    sin_half_phi = math.sin(delta_phi / 2.0)
    sin_half_lambda = math.sin(delta_lambda / 2.0)
    h = (
        sin_half_phi * sin_half_phi
        + math.cos(phi1) * math.cos(phi2) * sin_half_lambda * sin_half_lambda
    )
    h = max(0.0, min(1.0, h))
    return 2.0 * R * math.asin(math.sqrt(h))


def _point_in_polygon(x: float, y: float, poly: Sequence[tuple[float, float]]) -> bool:
    """Even-odd rule ray casting for point in polygon (x, y)."""
    inside = False
    n = len(poly)
    p1x, p1y = poly[0]
    for i in range(1, n + 1):
        p2x, p2y = poly[i % n]
        if min(p1y, p2y) < y <= max(p1y, p2y) and p1y != p2y:
            x_inters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
            if x <= x_inters:
                inside = not inside
        p1x, p1y = p2x, p2y
    return inside


def polygon_to_cells(
    ring: Sequence[tuple[float, float]],
    res: int,
    *,
    max_cells: int = 100_000,
) -> list[int]:
    """Return cell identifiers whose projected centres lie inside the ring."""
    if len(ring) > MAX_RING_VERTICES:
        raise GeoHexError("input_too_large", f"Ring exceeds maximum {MAX_RING_VERTICES} vertices")
    if len(ring) < 3:
        return []
    if not (0 <= res <= MAX_RESOLUTION):
        raise GeoHexError("invalid_resolution", f"Resolution must be 0..{MAX_RESOLUTION}")

    min_lng = min(pt[1] for pt in ring)
    max_lng = max(pt[1] for pt in ring)
    if max_lng - min_lng > 180.0:
        raise GeoHexError(
            "antimeridian_not_supported",
            "Ring spanning more than 180 degrees of longitude is not supported",
        )

    # Project vertices
    proj_ring: list[tuple[float, float]] = []
    for lat, lng in ring:
        proj_ring.append(project(lat, lng))

    min_x = min(pt[0] for pt in proj_ring)
    max_x = max(pt[0] for pt in proj_ring)
    min_y = min(pt[1] for pt in proj_ring)
    max_y = max(pt[1] for pt in proj_ring)

    # Find axial bounds over the 4 corners of bounding box
    corners = [
        (min_x, min_y),
        (min_x, max_y),
        (max_x, min_y),
        (max_x, max_y),
    ]
    corner_qrs = [_xy_to_qr(x, y, res) for x, y in corners]
    min_q = min(c[0] for c in corner_qrs) - 2
    max_q = max(c[0] for c in corner_qrs) + 2
    min_r = min(c[1] for c in corner_qrs) - 2
    max_r = max(c[1] for c in corner_qrs) + 2

    matched_cells: list[int] = []
    for q in range(min_q, max_q + 1):
        for r in range(min_r, max_r + 1):
            cx, cy = _qr_to_xy(q, r, res)
            if (
                min_x <= cx <= max_x
                and min_y <= cy <= max_y
                and _point_in_polygon(cx, cy, proj_ring)
            ):
                matched_cells.append(_encode_cell(res, q, r))
                if len(matched_cells) > max_cells:
                    raise GeoHexError(
                        "too_many_cells",
                        f"Exceeded max_cells limit ({max_cells})",
                    )

    return matched_cells
