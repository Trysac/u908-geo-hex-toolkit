# Aperture-7 Equal-Area Hexagonal Grid using Lambert Cylindrical Projection

Status: accepted (generated, pending human review)

## Context

Spatial analysis pipelines often require equal-area spatial binning without the binary dependency footprint of C extensions like H3 or GDAL, while maintaining deterministic hierarchical parent-child relationships.

## Decision

We implemented an aperture-7 discrete hexagonal grid system projected via the Lambert cylindrical equal-area projection using the authalic radius. Axial coordinates are packed into 64-bit integers with exact algebraic parent-child transitions.

## Consequences

Provides pure-Python equal-area bins with zero runtime dependencies. The grid does not wrap seamlessly at the antimeridian, rings spanning over 180 degrees of longitude are explicitly rejected with 'antimeridian_not_supported'.
