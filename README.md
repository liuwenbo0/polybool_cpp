# polybool_cpp

C++17 port of a polygon boolean operation library.

## Files

- `Polybool.hpp` and `Polybool.cpp`: polygon boolean implementation.
- `test.cpp`: standalone regression tests for polygon operations and shape output.

## Build And Test

```bash
g++ -std=c++17 -Wall -Wextra -pedantic test.cpp Polybool.cpp -o polybool_test
./polybool_test
```

Expected result:

```text
Pass: 33
Fail: 0
```

## Covered Operations

The test suite covers:

- union
- intersection
- difference
- reverse difference
- xor
- disjoint, overlapping, identical, edge-touching, and point-touching polygons
- contained polygons and hole-like region output
- multi-region polygons
- inverted polygons
- cubic curve output
- shape transforms
