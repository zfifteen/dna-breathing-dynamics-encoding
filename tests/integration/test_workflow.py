"""
Integration tests for end-to-end workflows.

Tests complete data processing pipelines combining multiple components.
"""

import pytest

from src.core.params import validate_dna_sequence, validate_temperature


@pytest.mark.integration
class TestBasicWorkflow:
    """Test basic analysis workflow."""

    def test_sequence_validation_workflow(
        self, long_dna_sequence: str, default_temperature: float
    ) -> None:
        """Test complete sequence validation workflow."""
        # Step 1: Validate sequence
        validated_seq = validate_dna_sequence(long_dna_sequence)
        assert len(validated_seq) == len(long_dna_sequence)
        assert validated_seq.isupper()

        # Step 2: Validate temperature
        validated_temp = validate_temperature(default_temperature)
        assert validated_temp == default_temperature

        # Step 3: Verify sequence properties
        assert all(base in "ACGTN" for base in validated_seq)

    def test_gc_content_workflow(self, gc_rich_sequence: str) -> None:
        """Test GC content analysis workflow."""
        validated_seq = validate_dna_sequence(gc_rich_sequence)

        # Calculate GC content
        g_count = validated_seq.count("G")
        c_count = validated_seq.count("C")
        gc_content = (g_count + c_count) / len(validated_seq)

        assert gc_content > 0.8  # GC-rich sequence


@pytest.mark.integration
@pytest.mark.smoke
def test_quick_integration_check(sample_dna_sequence: str) -> None:
    """Quick integration smoke test (<1s)."""
    # This should complete very quickly
    result = validate_dna_sequence(sample_dna_sequence)
    assert result is not None
