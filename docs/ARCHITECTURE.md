# Architecture

`geo_hex_toolkit` implements an aperture-7 hexagonal discrete global grid system (DGGS) using standard Python without external runtime dependencies.

## Projection and Lattice

- **Projection**: Lambert cylindrical equal-area projection:
  $$
  x = R \cdot \text{radians}(\text{lng}), \quad y = R \cdot \sin(\text{radians}(\text{lat}))
  $$
  using the authalic sphere radius $R = 6371007.180918475\text{ m}$.
- **Lattice**: Resolution 0 consists of pointy-top hexagons with edge length $e_0 = 1,000,000\text{ m}$. The axial basis vectors are:
  $$
  b_1 = (\sqrt{3} e_0, 0), \quad b_2 = (\frac{\sqrt{3}}{2} e_0, 1.5 e_0)
  $$
- **Hierarchy**: Each resolution step scales the grid by an aperture of 7 (area reduced by factor 7, edge length shrunk by $\sqrt{7}$, rotated by $\arctan(\sqrt{3}/5)$).
- **Parent / Child Links**: Resolution transition matrix $M = \begin{pmatrix} 2 & -1 \\ 1 & 3 \end{pmatrix}$ maps coarse axial unit vectors to fine axial coordinates $(2, 1)$ and $(-1, 3)$. Parent resolution rounding maps fine $(a, b)$ to cube rounding of $((3a + b)/7, (-a + 2b)/7)$.

## Cell Identifiers

Cell IDs are packed into 64-bit unsigned integers:

```
[63..60] Resolution (4 bits, 0..15)
[59..30] Axial q + 2^29 (30 bits)
[29..0]  Axial r + 2^29 (30 bits)
```

Coordinates outside $[-2^{29}, 2^{29})$ raise `GeoHexError("invalid_cell")`.
