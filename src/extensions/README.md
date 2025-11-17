# C/C++ Extensions

High-performance extensions for DNA breathing dynamics analysis with Apple Silicon optimization.

## Overview

This directory contains C/C++ extensions for performance-critical operations:
- DNA sequence analysis and validation
- High-precision numerical computations
- SIMD-optimized algorithms
- Apple Silicon AMX instruction utilization

## Building Extensions

To build the C/C++ extensions:

```bash
make build
```

Or directly:

```bash
python setup.py build_ext --inplace
```

## Apple Silicon Optimization

Extensions are optimized for Apple Silicon (M1/M2/M3/M4) with:
- **AMX instructions** for matrix operations
- **SIMD vectorization** (NEON)
- **Cache-line alignment** for memory operations
- **Metal compute shader support** (future)

## File Structure

- `common.h` - Shared definitions, platform detection, utilities
- `fast_compute_example.c` - Example extension (template)
- Additional extensions can be added as needed

## Adding New Extensions

1. Create new `.c` file with Python C API bindings
2. Include `common.h` for shared utilities
3. Add extension definition to `setup.py`
4. Rebuild with `make build`

## Dependencies

Required system libraries (install via `make deps`):
- **MPFR** - High-precision floating point
- **GMP** - Arbitrary precision arithmetic
- **libomp** - OpenMP support (optional)

## Platform Detection

Extensions automatically detect Apple Silicon:

```c
#ifdef APPLE_SILICON
    // Apple Silicon-specific code
#endif
```

## Performance Considerations

- Use `INLINE` for frequently called functions
- Apply `ALIGNED(64)` for cache-line alignment
- Use `LIKELY`/`UNLIKELY` for branch prediction hints
- Leverage `RESTRICT` for pointer aliasing
