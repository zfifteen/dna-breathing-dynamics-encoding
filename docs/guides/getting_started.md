# Getting Started

This guide will help you get up and running with the DNA Breathing Dynamics framework.

## Prerequisites

### System Requirements

- **macOS** (Apple Silicon preferred) or **Linux**
- **Python** 3.10 or higher
- **Homebrew** (macOS) or appropriate package manager
- 4GB+ RAM (8GB+ recommended for large-scale analysis)

### Required System Libraries

**macOS (via Homebrew):**
```bash
brew install mpfr gmp libomp
```

**Linux (Ubuntu/Debian):**
```bash
sudo apt-get install libmpfr-dev libgmp-dev libomp-dev
```

## Installation

### 1. Clone the Repository

```bash
git clone <repository-url>
cd dna_breathing_dynamics_encoding
```

### 2. Install Dependencies

Using the Makefile (recommended):
```bash
# Install system dependencies (macOS only)
make deps

# Create virtual environment and install Python packages
make install
```

Manual installation:
```bash
# Create virtual environment
python3 -m venv .venv
source .venv/bin/activate

# Install Python dependencies
pip install -r requirements.txt
pip install -e .
```

### 3. Build C/C++ Extensions

```bash
make build
```

Or manually:
```bash
python setup.py build_ext --inplace
```

### 4. Verify Installation

```bash
# Run smoke tests (<5s)
make smoke
```

Expected output:
```
==================== test session starts ====================
collected 3 items

tests/unit/test_params.py::test_parameter_defaults PASSED
tests/integration/test_workflow.py::test_quick_integration_check PASSED

===================== 2 passed in 0.5s =====================
```

## Basic Usage

### Validating DNA Sequences

```python
from src.core.params import validate_dna_sequence

# Valid DNA sequence
sequence = "ATCGATCGATCG"
validated = validate_dna_sequence(sequence)
print(validated)  # ATCGATCGATCG

# Automatic uppercase conversion
sequence_lower = "atcg"
validated = validate_dna_sequence(sequence_lower)
print(validated)  # ATCG

# RNA sequence (with U instead of T)
rna_sequence = "AUCGAUCG"
validated_rna = validate_dna_sequence(rna_sequence, allow_rna=True)
print(validated_rna)  # AUCGAUCG
```

### Using Default Parameters

```python
from src.core.params import (
    TEMPERATURE_DEFAULT,
    KAPPA_GEO_DEFAULT,
    validate_temperature,
)

# Default physiological temperature
temp = TEMPERATURE_DEFAULT
print(f"Temperature: {temp}K ({temp - 273.15}°C)")
# Temperature: 310.15K (37.0°C)

# Validate custom temperature
custom_temp = validate_temperature(298.15)  # 25°C
```

### High-Precision Calculations

```python
from mpmath import mp
from src.core.params import DEFAULT_MPMATH_DPS

# Set precision (already configured automatically)
mp.dps = DEFAULT_MPMATH_DPS

# High-precision calculation example
from mpmath import mpf, exp

result = exp(mpf("100"))
print(f"e^100 = {result}")
```

## Running Tests

### Quick Validation

```bash
# Smoke tests (<5s)
make smoke
```

### Comprehensive Testing

```bash
# All tests
make test

# Specific test categories
make unit           # Unit tests only
make integration    # Integration tests
make performance    # Performance benchmarks
```

### With Coverage

```bash
pytest tests/ -v --cov=src --cov-report=html
open htmlcov/index.html  # View coverage report
```

## Development Workflow

### 1. Activate Virtual Environment

```bash
source .venv/bin/activate
```

### 2. Make Changes

Edit files in `src/`, `tests/`, etc.

### 3. Format and Lint

```bash
make format    # Auto-format with Black and isort
make lint      # Check with Flake8
make typecheck # Type check with MyPy
```

### 4. Run Tests

```bash
make test
```

### 5. Build and Test Extensions (if modified)

```bash
make clean
make build
make test
```

## Project Commands

### Build Commands

```bash
make deps      # Install system dependencies (macOS)
make venv      # Create virtual environment
make install   # Install Python dependencies
make build     # Build C/C++ extensions
make clean     # Remove build artifacts
```

### Testing Commands

```bash
make test          # Run all tests
make smoke         # Quick smoke tests (<5s)
make unit          # Unit tests only
make integration   # Integration tests
make performance   # Performance benchmarks
make validation    # Scientific validation tests
make bench         # Run benchmarks with autosave
```

### Code Quality Commands

```bash
make format    # Format code (Black + isort)
make lint      # Lint with Flake8
make typecheck # Type check with MyPy
make check     # All quality checks
```

### Comprehensive Commands

```bash
make all       # Full build and test cycle
make ci        # CI/CD pipeline (no deps install)
make dev       # Set up development environment
```

## Next Steps

- Read the [API Documentation](../api/) for detailed module information
- Explore [Examples](../../examples/) for common use cases
- Review [Validation Guide](./validation.md) for scientific rigor standards
- Check [Performance Optimization](./performance.md) for Apple Silicon tips

## Troubleshooting

### Build Errors

**Issue**: C extension compilation fails
```bash
# Verify system libraries are installed
brew list | grep -E "(mpfr|gmp|openmp)"  # macOS
ldconfig -p | grep -E "(mpfr|gmp)"       # Linux
```

**Issue**: Wrong Homebrew prefix
```bash
# Check your Homebrew location
which brew
echo $PATH | tr ':' '\n' | grep homebrew
```

### Import Errors

**Issue**: Cannot import src modules
```bash
# Ensure package is installed in editable mode
pip install -e .
```

### Test Failures

**Issue**: Precision-related test failures
```python
# Verify mpmath precision settings
import mpmath
print(mpmath.mp.dps)  # Should be 50
```

## Getting Help

- Check documentation in `docs/`
- Review existing tests in `tests/` for examples
- Open an issue on GitHub
- Consult CLAUDE.md for project-specific guidance
