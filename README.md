# DNA Breathing Dynamics and Encoding

A hybrid bioinformatics and mathematical modeling framework for analyzing DNA breathing dynamics with high-precision numerical computation, optimized for Apple Silicon.

## Overview

This framework combines biological sequence analysis with mathematical modeling to study DNA breathing dynamics - the transient opening and closing of base pairs in the double helix. It provides:

- **Bioinformatics Analysis**: DNA/RNA sequence processing with strict validation
- **Mathematical Modeling**: High-precision numerical computations (mpmath, MPFR)
- **Apple Silicon Optimization**: AMX instructions and SIMD optimizations for M1/M2/M3 Macs
- **Scientific Rigor**: Statistical validation, bootstrap confidence intervals, empirical testing

## Quick Start

### Installation

```bash
# Install system dependencies (macOS)
make deps

# Create virtual environment and install Python dependencies
make install

# Build C/C++ extensions
make build
```

### Running Tests

```bash
# Quick smoke test (<5s)
make smoke

# Full test suite
make test

# Performance benchmarks
make bench
```

### Basic Usage

```python
from src.core.params import validate_dna_sequence, TEMPERATURE_DEFAULT

# Validate DNA sequence
sequence = "ATCGATCGATCG"
validated = validate_dna_sequence(sequence)

# Use default parameters
temp = TEMPERATURE_DEFAULT  # 310.15 K (37°C)
```

## Project Structure

```
dna_breathing_dynamics_encoding/
├── src/
│   ├── core/          # Core framework components
│   ├── bio/           # Bioinformatics modules
│   ├── math/          # Mathematical modeling
│   └── extensions/    # C/C++ high-performance extensions
├── tests/
│   ├── unit/          # Unit tests
│   ├── integration/   # Integration tests
│   └── performance/   # Performance benchmarks
├── proof_pack/        # Scientific validation framework
├── examples/          # Usage examples
├── docs/              # Documentation
└── scripts/           # Utility scripts
```

## Features

### Bioinformatics

- DNA/RNA sequence validation with nucleotide checking
- Graph-theoretic validation for dinucleotide-preserving shuffles (Eulerian path analysis)
- GC content analysis
- Sequence property calculations
- Strict validation following human genome standards (GRCh38/hg38)

### Mathematical Framework

- High-precision arithmetic (mpmath dps=50, MPFR 256-bit)
- Geodesic mapping with κ_geo parameter
- Z_5D framework integration (κ_star calibration)
- Numerical stability monitoring

### Apple Silicon Optimization

- AMX instruction utilization for matrix operations
- SIMD vectorization (NEON)
- Cache-line aligned memory operations
- Platform-specific compiler optimizations

### Scientific Validation

- Statistical hypothesis testing (t-test, Wilcoxon, KS, χ²)
- Bootstrap confidence intervals (1000+ resamples)
- Numerical accuracy validation
- Convergence monitoring for iterative algorithms

## Development

### Setup Development Environment

```bash
make dev
source .venv/bin/activate
```

### Code Quality

```bash
# Format code
make format

# Lint code
make lint

# Type check
make typecheck

# Run all checks
make check
```

### Building Documentation

```bash
make docs
```

## Testing Strategy

### Test Categories

- **Unit tests** (`tests/unit/`): Component-level testing
- **Integration tests** (`tests/integration/`): End-to-end workflows
- **Performance tests** (`tests/performance/`): Benchmarks and regression testing
- **Validation tests** (`proof_pack/`): Scientific validation and statistical rigor

### Running Specific Tests

```bash
# Unit tests only
make unit

# Integration tests
make integration

# Performance benchmarks
make performance

# Scientific validation
make validation
```

## Scientific Rigor

This framework follows strict scientific computing standards:

- **Exact version pinning** for reproducibility
- **Bootstrap confidence intervals** for all statistical claims
- **Empirical validation** over hypothetical extrapolations
- **Cross-platform compatibility** with documented environments
- **Pre-registered test endpoints** to prevent data leakage

## Platform Support

- **Python**: 3.10, 3.11, 3.12
- **Operating Systems**: macOS (Apple Silicon preferred), Linux
- **Required Libraries**: MPFR, GMP, libomp (via Homebrew)

## Performance Considerations

### Apple Silicon Optimizations

When running on M1/M2/M3/M4 Macs:
- Automatic detection of Apple Silicon architecture
- AMX instruction utilization for matrix operations
- Optimized Homebrew library paths (`/opt/homebrew`)
- High memory capacity utilization for large-scale computations

### Precision Settings

- Default mpmath precision: 50 decimal places
- Default MPFR precision: 256 bits
- Configurable via `src/core/params.py`

## Contributing

### Code Style

- **Black** formatting (88 character line limit)
- **isort** import organization
- **Flake8** linting compliance
- **MyPy** type annotations

### Adding New Features

1. Create feature branch from `main`
2. Implement with tests (unit + integration)
3. Add documentation
4. Run full validation: `make all`
5. Submit pull request with CI passing

## License

MIT License - See LICENSE file for details

## References

- Mathematical framework patterns from unified-framework
- Bioinformatics standards from CRISPR/DNA analysis projects
- Apple Silicon optimization techniques from z-amx research

## Acknowledgments

Built following best practices from the VelocityWorks project portfolio, including rigorous scientific validation, Apple Silicon optimization, and comprehensive testing standards.
