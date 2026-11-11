# Project record

## Objective

Pure Python hierarchical hexagon grid on an equal-area projection with exact integer parent and child links, great-circle distances and polygon cover

## Scope

- In scope: the public API documented in [README.md](../README.md), typed and tested.
- Out of scope: anything that needs runtime dependencies or I/O without an accepted ADR.

## Timeline

| Date | Event |
|---|---|
| 2026 | Repository created from the factory utility template. |

## Effort metrics

Measured values only; a value that has not been measured stays "pending measurement".

| Metric | Value |
|---|---|
| Tokens (input / output) | 2233715 / 25352 |
| Hours (human) | pending measurement |
| Cost (USD) | 0.8142 |

## Factory run u908 (measured from API usage)

| Phase | Model | Batch | Input | Output | Cache write | Cache read | Cost (USD) |
|---|---|---|---|---|---|---|---|
| generation | google/gemini-3.8-flash | yes | 5934 | 10402 | 0 | 0 | 0.0435 |
| fix-1 | google/gemini-3.8-flash | no | 200965 | 7464 | 0 | 404241 | 0.2090 |
| fix-2 | google/gemini-3.8-flash | no | 249504 | 3061 | 0 | 371896 | 0.2265 |
| fix-3 | google/gemini-3.8-flash | no | 151545 | 2501 | 0 | 363733 | 0.1503 |
| fix-4 | google/gemini-3.8-flash | no | 204855 | 1298 | 0 | 94512 | 0.1656 |
| audit | google/gemini-3.1-flash-lite | no | 60775 | 626 | 0 | 125755 | 0.0193 |

Total: $0.8142. Prices: vertex line, as measured by its adapter.
