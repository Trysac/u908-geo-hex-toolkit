# geo-hex-toolkit

[![PyPI version](https://img.shields.io/pypi/v/trysac-geo-hex-toolkit)](https://pypi.org/project/trysac-geo-hex-toolkit/)
[![License: MIT](https://img.shields.io/pypi/l/trysac-geo-hex-toolkit)](LICENSE)
[![CI](https://github.com/Trysac/u908-geo-hex-toolkit/actions/workflows/ci.yml/badge.svg)](https://github.com/Trysac/u908-geo-hex-toolkit/actions/workflows/ci.yml)

Pure Python hierarchical hexagon grid on an equal-area projection with exact integer parent and child links, great-circle distances and polygon cover

## Problem

Describe in two or three sentences the concrete problem this utility solves, who hits it, and what
they do today without it. Replace this paragraph when the first real feature lands.

## Install

```sh
pip install trysac-geo-hex-toolkit
# or
uv add trysac-geo-hex-toolkit
```

Requires Python 3.10 or newer. No runtime dependencies.

## Usage

```python
import geo_hex_toolkit

# normalize() trims surrounding whitespace and lower-cases the text.
print(geo_hex_toolkit.normalize("  Hello World "))  # -> "hello world"

# Blank input is rejected explicitly instead of returning an empty string.
try:
    geo_hex_toolkit.normalize("   ")
except ValueError as exc:
    print(f"rejected: {exc}")
```

## Documentation

- [Architecture](docs/ARCHITECTURE.md) and [decisions](docs/decisions/)
- [Testing](docs/TESTING.md) and [security](docs/SECURITY.md)
- [Project record](docs/PROJECT.md) and [stack](docs/STACK.md)
- [Changelog](CHANGELOG.md) and [contributing guide](CONTRIBUTING.md)

## License

[MIT](LICENSE) (c) 2026 Miguel Gomez
