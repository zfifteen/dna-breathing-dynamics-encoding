# API Overview

Comprehensive overview of the DNA Breathing Dynamics framework API.

## Module Structure

### Core Modules (`src/core/`)

#### `src.core.params`

Centralized parameter management and validation.

**Constants:**
- `DEFAULT_MPMATH_DPS` (int): Default mpmath decimal precision (50)
- `DEFAULT_MPFR_PRECISION` (int): Default MPFR bit precision (256)
- `DNA_NUCLEOTIDES` (frozenset): Valid DNA nucleotide characters
- `RNA_NUCLEOTIDES` (frozenset): Valid RNA nucleotide characters
- `TEMPERATURE_DEFAULT` (float): Default temperature in Kelvin (310.15 K)
- `KAPPA_GEO_DEFAULT` (float): Default geodesic exponent (0.3)
- `KAPPA_STAR_DEFAULT` (float): Default Z_5D calibration (0.04449)

**Functions:**

```python
validate_dna_sequence(
    sequence: str,
    allow_rna: bool = False,
    allow_lowercase: bool = True
) -> str
```
Validate DNA/RNA sequence composition.

```python
validate_temperature(temp: float) -> float
```
Validate temperature parameter (273.15 - 373.15 K).

```python
validate_kappa_geo(kappa_geo: float) -> float
```
Validate geodesic exponent parameter (0.0 - 1.0).

```python
validate_kappa_star(kappa_star: float) -> float
```
Validate Z_5D calibration parameter (0.0 - 0.1).

```python
validate_positive(value: float, name: str = "value") -> float
```
Validate positive numeric parameter.

```python
validate_probability(prob: float, name: str = "probability") -> float
```
Validate probability in [0, 1].

### Bioinformatics Modules (`src/bio/`)

**Status**: Ready for implementation

Planned modules:
- `sequence.py`: DNA/RNA sequence analysis
- `structure.py`: Secondary structure prediction
- `breathing.py`: Breathing dynamics characterization
- `thermodynamics.py`: Base pair stability calculations

### Mathematical Modules (`src/math/`)

**Status**: Ready for implementation

Planned modules:
- `precision.py`: High-precision arithmetic utilities
- `geodesic.py`: Geodesic mapping implementations
- `statistics.py`: Statistical analysis functions
- `optimization.py`: Numerical optimization algorithms

### Extensions (`src/extensions/`)

High-performance C/C++ extensions with Apple Silicon optimization.

**Common Utilities** (`common.h`):
- Platform detection (APPLE_SILICON)
- Optimization macros (INLINE, RESTRICT, LIKELY, UNLIKELY)
- Nucleotide encoding (encode_nucleotide, is_gc, is_at)

**Example Extensions** (`fast_compute_example.c`):
- `fast_gc_content(sequence: str) -> float`
- `fast_validate_sequence(sequence: str) -> bool`

## Validation Framework (`proof_pack/`)

### `proof_pack.validation_framework`

Scientific validation with statistical rigor.

#### Classes

**ValidationResult**

Container for validation test results.

Attributes:
- `test_name` (str): Name of the test
- `passed` (bool): Whether test passed
- `p_value` (float | None): Statistical p-value
- `statistic` (float | None): Test statistic
- `confidence_interval` (tuple | None): Bootstrap CI
- `effect_size` (float | None): Effect size (Cohen's d)
- `metadata` (dict): Additional test metadata

Methods:
- `summary() -> str`: Human-readable summary

**ValidationSuite**

Orchestrate multiple validation tests.

Methods:
- `add_result(result: ValidationResult)`: Add test result
- `run_all() -> bool`: Check if all tests passed
- `summary() -> str`: Generate summary report
- `detailed_report() -> str`: Generate detailed report

#### Functions

```python
bootstrap_confidence_interval(
    data: np.ndarray,
    statistic_func: Callable,
    n_resamples: int = 1000,
    confidence_level: float = 0.95,
    random_seed: int | None = None
) -> tuple[float, float]
```
Calculate bootstrap confidence interval.

```python
validate_statistical_hypothesis(
    observed_data: np.ndarray,
    expected_data: np.ndarray | None = None,
    test_type: str = "t-test",
    alpha: float = 0.05,
    alternative: str = "two-sided"
) -> ValidationResult
```
Perform statistical hypothesis testing.

Supported test types:
- `"t-test"`: Student's t-test
- `"wilcoxon"`: Wilcoxon signed-rank test
- `"ks"`: Kolmogorov-Smirnov test
- `"chi2"`: Chi-squared test

```python
validate_numerical_accuracy(
    computed_values: np.ndarray,
    reference_values: np.ndarray,
    rtol: float = 1e-9,
    atol: float = 1e-12
) -> ValidationResult
```
Validate numerical accuracy against reference values.

```python
validate_convergence(
    sequence: np.ndarray,
    tolerance: float = 1e-10,
    window: int = 10
) -> ValidationResult
```
Validate convergence of iterative algorithm.

## Usage Examples

### Parameter Validation

```python
from src.core.params import (
    validate_dna_sequence,
    validate_temperature,
    TEMPERATURE_DEFAULT
)

# Validate sequence
seq = validate_dna_sequence("ATCG")

# Validate temperature
temp = validate_temperature(310.15)

# Use defaults
default_temp = TEMPERATURE_DEFAULT
```

### Statistical Validation

```python
from proof_pack import (
    ValidationSuite,
    validate_statistical_hypothesis,
    bootstrap_confidence_interval
)
import numpy as np

# Create validation suite
suite = ValidationSuite("My Analysis")

# Statistical test
observed = np.random.normal(1.0, 0.5, 100)
expected = np.random.normal(0.0, 0.5, 100)

result = validate_statistical_hypothesis(
    observed, expected, test_type="t-test"
)
suite.add_result(result)

# Bootstrap CI
data = np.random.exponential(2.0, 200)
ci = bootstrap_confidence_interval(
    data, np.mean, n_resamples=1000
)

# Summary
print(suite.summary())
```

### High-Precision Computation

```python
from mpmath import mp, mpf, exp
from src.core.params import DEFAULT_MPMATH_DPS

# Configure precision
mp.dps = DEFAULT_MPMATH_DPS

# High-precision calculation
x = mpf("0.1")
result = exp(x)
print(f"e^0.1 = {result}")
```

## Type Annotations

All public APIs include comprehensive type annotations compatible with MyPy.

Example:
```python
def validate_dna_sequence(
    sequence: str,
    allow_rna: bool = False,
    allow_lowercase: bool = True
) -> str:
    ...
```

## Error Handling

All validation functions raise `ValueError` with descriptive messages:

```python
try:
    validate_dna_sequence("ATCGXYZ")
except ValueError as e:
    print(e)  # Invalid nucleotide characters: ['X', 'Y', 'Z']...
```

## Extension API (C/C++)

### Error Codes

```c
typedef enum {
    DNA_SUCCESS = 0,
    DNA_ERROR_INVALID_SEQUENCE = -1,
    DNA_ERROR_INVALID_PARAMETER = -2,
    DNA_ERROR_MEMORY_ALLOCATION = -3,
    DNA_ERROR_COMPUTATION = -4,
} dna_error_t;
```

### Nucleotide Encoding

```c
typedef enum {
    NUC_A = 0, NUC_C = 1, NUC_G = 2,
    NUC_T = 3, NUC_N = 4, NUC_INVALID = -1
} nucleotide_t;

nucleotide_t encode_nucleotide(char c);
bool is_valid_dna_nucleotide(char c);
bool is_gc(nucleotide_t nuc);
bool is_at(nucleotide_t nuc);
```

## Platform-Specific Features

### Apple Silicon Detection

```python
# Python
import platform
is_apple_silicon = (
    platform.system() == "Darwin" and
    platform.machine() == "arm64"
)
```

```c
// C/C++
#ifdef APPLE_SILICON
    // Apple Silicon specific code
#endif
```

## Performance Considerations

- **Sequence validation**: O(n) where n = sequence length
- **Statistical tests**: O(n) to O(n²) depending on test type
- **Bootstrap CI**: O(n × k) where k = n_resamples
- **High-precision arithmetic**: Varies by operation and precision

## Thread Safety

- Core parameter validation: Thread-safe
- mpmath operations: Not thread-safe (use process-level parallelism)
- C extensions: Thread-safe where documented

## Future API Extensions

Planned additions:
- DNA breathing dynamics simulation
- Free energy calculations
- Machine learning integration
- GPU acceleration (Metal compute shaders)
