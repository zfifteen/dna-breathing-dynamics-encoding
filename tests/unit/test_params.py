"""
Unit tests for core parameter management and validation.

Tests parameter validation functions, bounds checking, and error handling.
"""

import pytest

from src.core.params import (
    DNA_NUCLEOTIDES,
    KAPPA_GEO_DEFAULT,
    KAPPA_STAR_DEFAULT,
    RNA_NUCLEOTIDES,
    TEMPERATURE_DEFAULT,
    validate_dna_sequence,
    validate_kappa_geo,
    validate_kappa_star,
    validate_positive,
    validate_probability,
    validate_temperature,
)


@pytest.mark.unit
class TestDNASequenceValidation:
    """Test DNA/RNA sequence validation."""

    def test_valid_dna_sequence(self, sample_dna_sequence: str) -> None:
        """Test validation of valid DNA sequence."""
        result = validate_dna_sequence(sample_dna_sequence)
        assert result == sample_dna_sequence.upper()

    def test_valid_rna_sequence(self, sample_rna_sequence: str) -> None:
        """Test validation of valid RNA sequence."""
        result = validate_dna_sequence(sample_rna_sequence, allow_rna=True)
        assert result == sample_rna_sequence.upper()

    def test_lowercase_conversion(self) -> None:
        """Test automatic uppercase conversion."""
        seq = "atcg"
        result = validate_dna_sequence(seq)
        assert result == "ATCG"

    def test_invalid_characters(self) -> None:
        """Test rejection of invalid nucleotide characters."""
        with pytest.raises(ValueError, match="Invalid nucleotide characters"):
            validate_dna_sequence("ATCGXYZ")

    def test_empty_sequence(self) -> None:
        """Test rejection of empty sequences."""
        with pytest.raises(ValueError, match="Sequence cannot be empty"):
            validate_dna_sequence("")

    def test_rna_uracil_without_flag(self) -> None:
        """Test that U is rejected without allow_rna flag."""
        with pytest.raises(ValueError):
            validate_dna_sequence("AUCG", allow_rna=False)


@pytest.mark.unit
class TestTemperatureValidation:
    """Test temperature parameter validation."""

    def test_valid_temperature(self, default_temperature: float) -> None:
        """Test validation of valid temperature."""
        result = validate_temperature(default_temperature)
        assert result == default_temperature

    def test_temperature_too_low(self) -> None:
        """Test rejection of temperature below minimum."""
        with pytest.raises(ValueError, match="outside valid range"):
            validate_temperature(200.0)  # Below 273.15K

    def test_temperature_too_high(self) -> None:
        """Test rejection of temperature above maximum."""
        with pytest.raises(ValueError, match="outside valid range"):
            validate_temperature(400.0)  # Above 373.15K

    def test_boundary_temperatures(self) -> None:
        """Test boundary temperature values."""
        assert validate_temperature(273.15) == 273.15
        assert validate_temperature(373.15) == 373.15


@pytest.mark.unit
class TestKappaValidation:
    """Test kappa parameter validation."""

    def test_valid_kappa_geo(self) -> None:
        """Test validation of valid kappa_geo."""
        result = validate_kappa_geo(KAPPA_GEO_DEFAULT)
        assert result == KAPPA_GEO_DEFAULT

    def test_kappa_geo_out_of_range(self) -> None:
        """Test rejection of kappa_geo outside valid range."""
        with pytest.raises(ValueError):
            validate_kappa_geo(-0.1)
        with pytest.raises(ValueError):
            validate_kappa_geo(1.5)

    def test_valid_kappa_star(self) -> None:
        """Test validation of valid kappa_star."""
        result = validate_kappa_star(KAPPA_STAR_DEFAULT)
        assert result == KAPPA_STAR_DEFAULT

    def test_kappa_star_out_of_range(self) -> None:
        """Test rejection of kappa_star outside valid range."""
        with pytest.raises(ValueError):
            validate_kappa_star(-0.01)
        with pytest.raises(ValueError):
            validate_kappa_star(0.15)


@pytest.mark.unit
class TestGenericValidation:
    """Test generic validation functions."""

    def test_validate_positive(self) -> None:
        """Test positive value validation."""
        assert validate_positive(1.0) == 1.0
        assert validate_positive(0.001) == 0.001

        with pytest.raises(ValueError, match="must be positive"):
            validate_positive(0.0)

        with pytest.raises(ValueError, match="must be positive"):
            validate_positive(-1.0)

    def test_validate_probability(self) -> None:
        """Test probability validation."""
        assert validate_probability(0.0) == 0.0
        assert validate_probability(0.5) == 0.5
        assert validate_probability(1.0) == 1.0

        with pytest.raises(ValueError, match="must be in"):
            validate_probability(-0.1)

        with pytest.raises(ValueError, match="must be in"):
            validate_probability(1.5)


@pytest.mark.unit
@pytest.mark.smoke
def test_parameter_defaults() -> None:
    """Quick smoke test for parameter default values."""
    assert KAPPA_GEO_DEFAULT == 0.3
    assert KAPPA_STAR_DEFAULT == 0.04449
    assert TEMPERATURE_DEFAULT == 310.15
    assert "A" in DNA_NUCLEOTIDES
    assert "U" in RNA_NUCLEOTIDES
