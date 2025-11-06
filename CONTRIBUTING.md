# Contributing to DNA Breathing Dynamics Framework

Thank you for your interest in contributing! This document provides guidelines and best practices for contributing to the project.

## Development Setup

1. **Fork and clone the repository**
   ```bash
   git clone <your-fork-url>
   cd dna_breathing_dynamics_encoding
   ```

2. **Set up development environment**
   ```bash
   make deps    # Install system dependencies (macOS)
   make install # Install Python dependencies
   make build   # Build C/C++ extensions
   ```

3. **Verify setup**
   ```bash
   make smoke   # Run quick tests
   ```

## Development Workflow

### 1. Create a Feature Branch

```bash
git checkout -b feature/your-feature-name
```

Use descriptive branch names:
- `feature/add-breathing-simulation`
- `fix/sequence-validation-bug`
- `docs/improve-api-documentation`
- `perf/optimize-gc-calculation`

### 2. Make Your Changes

Follow these guidelines:

#### Code Style

- **Python**: Follow PEP 8 with Black formatting (88 char line length)
- **C/C++**: 4-space indentation, descriptive variable names
- **Documentation**: Google-style docstrings

Run code quality checks:
```bash
make format    # Auto-format code
make lint      # Check linting
make typecheck # Type checking
```

#### Type Annotations

Always include type annotations for public APIs:
```python
def validate_sequence(sequence: str, allow_rna: bool = False) -> str:
    """Validate DNA/RNA sequence."""
    ...
```

#### Documentation

- Add docstrings to all public functions/classes
- Update relevant documentation in `docs/`
- Include usage examples where appropriate

### 3. Write Tests

Every new feature or bug fix must include tests:

```bash
# Create tests in appropriate directory
tests/unit/test_your_feature.py       # Unit tests
tests/integration/test_workflow.py     # Integration tests
tests/performance/test_benchmarks.py   # Performance tests
```

Test guidelines:
- Aim for >90% code coverage
- Include edge cases and error conditions
- Add smoke tests for critical paths
- Use descriptive test names

Example:
```python
@pytest.mark.unit
def test_sequence_validation_with_invalid_characters():
    """Test that invalid nucleotides are rejected."""
    with pytest.raises(ValueError, match="Invalid nucleotide"):
        validate_dna_sequence("ATCGXYZ")
```

### 4. Run Tests

```bash
make test      # All tests
make unit      # Unit tests only
make smoke     # Quick validation
```

Ensure all tests pass before submitting PR.

### 5. Commit Your Changes

Follow conventional commit format:

```
<type>(<scope>): <subject>

<body>

<footer>
```

Types:
- `feat`: New feature
- `fix`: Bug fix
- `docs`: Documentation changes
- `perf`: Performance improvements
- `refactor`: Code refactoring
- `test`: Test additions/changes
- `ci`: CI/CD changes

Examples:
```bash
git commit -m "feat(bio): Add GC content calculation function

Implement high-performance GC content analysis with SIMD
optimization for Apple Silicon.

Closes #42"
```

### 6. Push and Create Pull Request

```bash
git push origin feature/your-feature-name
```

Create a pull request with:
- Clear title and description
- Reference to related issues
- Screenshots/examples if applicable
- Confirmation that tests pass

## Code Review Process

1. **Automated Checks**: CI/CD runs automatically
   - All tests must pass
   - Code coverage must not decrease
   - Linting must pass

2. **Manual Review**: Maintainers will review for:
   - Code quality and style
   - Test coverage
   - Documentation completeness
   - Scientific rigor (if applicable)

3. **Iterations**: Address review feedback promptly

4. **Merge**: Once approved, maintainers will merge

## Scientific Contributions

For contributions involving scientific claims or algorithms:

### Requirements

1. **Empirical Validation**
   - Include validation tests in `proof_pack/`
   - Use bootstrap confidence intervals (n≥1000)
   - Document statistical methodology

2. **Reproducibility**
   - Set random seeds explicitly
   - Document all parameters
   - Include environment specifications

3. **Clear Labeling**
   - Distinguish validated results from hypotheses
   - Document assumptions clearly
   - Include references to literature

Example:
```python
from proof_pack import ValidationSuite, validate_statistical_hypothesis

suite = ValidationSuite("New Algorithm Validation")

# Empirical validation with statistical testing
result = validate_statistical_hypothesis(
    observed_data, expected_data,
    test_type="t-test", alpha=0.05
)
suite.add_result(result)

# Document in docstring
def new_algorithm(data: np.ndarray) -> float:
    """
    Novel algorithm for X analysis.

    Empirical Validation:
    - Validated on N=10^6 samples
    - Statistical significance: p < 0.001
    - Effect size: Cohen's d = 1.23
    - CI: [0.95, 1.01] (95% bootstrap, n=1000)

    Args:
        data: Input data array

    Returns:
        Calculated metric

    References:
        Smith et al. (2024). Journal of X. DOI:...
    """
    ...
```

## Performance Contributions

For performance optimizations:

1. **Benchmarks Required**
   ```python
   @pytest.mark.benchmark
   def test_optimized_function(benchmark):
       result = benchmark(optimized_function, data)
   ```

2. **Document Improvements**
   - Baseline performance
   - Optimized performance
   - Speedup factor
   - Platform tested (especially Apple Silicon)

3. **Maintain Correctness**
   - Numerical accuracy tests
   - Cross-platform compatibility

## Documentation Contributions

Documentation improvements are always welcome:

1. **API Documentation**: Update docstrings and `docs/api/`
2. **Guides**: Add tutorials in `docs/guides/`
3. **Examples**: Create demos in `examples/`
4. **README**: Improve main documentation

## Bug Reports

When reporting bugs, include:

1. **Description**: What happened vs. what you expected
2. **Reproduction**: Minimal code to reproduce
3. **Environment**:
   ```python
   import platform
   print(platform.platform())
   print(platform.python_version())
   ```
4. **Error Messages**: Full stack traces

## Feature Requests

When requesting features, include:

1. **Use Case**: What problem does this solve?
2. **Proposed Solution**: How might it work?
3. **Alternatives**: What alternatives have you considered?
4. **Impact**: Who benefits and how?

## Questions?

- Check existing documentation in `docs/`
- Search existing issues and PRs
- Review CLAUDE.md for project patterns
- Open a discussion issue

## Code of Conduct

- Be respectful and inclusive
- Focus on constructive feedback
- Assume good intentions
- Help newcomers

Thank you for contributing!
