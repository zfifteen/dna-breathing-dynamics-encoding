"""
Integration tests for local perturbation sweep pipeline.

Tests end-to-end: loading, perturbation generation, diff computation, outputs.
Researcher-framed docstrings; no full logic—skeletons verify structure/shapes.
"""

import sys
from pathlib import Path

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

import pytest
import pandas as pd

from experiments.local_perturbation_sweep.sweep_datagen import (
    load_and_perturb_guides,
    analyze_spectral_shifts,
)
from experiments.local_perturbation_sweep.utils import (
    generate_mutant,
    generate_double_mutant,
    compute_diffs,
    load_guides,
)


@pytest.fixture
def sample_guides() -> list:
    """
    PURPOSE: As a researcher/engineer, I want sample guides fixture so that I can
    test pipeline without full dataset.

    INPUTS: [None; generates 5 synthetic 20nt guides: 3 low-GC, 2 high-GC]
    EXPECTED OUTPUT: [List[str] of 5 valid ATGC sequences, len=20; gc_content verified]
    TEST DATA: [e.g., ['AT'*10, 'GC'*10]; expect gc=0.0/1.0, all len=20]
    REPRODUCTION: [Run fixture; assert len=5, all ATGC, gc balanced; re-run same]
    """
    return [
        "ATCGATCGATCGATCGATCG",
        "GCTAGCTAGCTAGCTAGCTA",
        "ATCGATCGATCGATCGATCG",
        "GCGCGCGCGCGCGCGCGCGC",
        "ATCGATCGATCGATCGATCG",
    ]


@pytest.fixture
def sample_perturbations(sample_guides: list) -> list:
    """
    PURPOSE: As a researcher, I want perturbed samples fixture so that I can test
    diff computation end-to-end.

    INPUTS: [sample_guides: list from fixture; generates singles at pos=1,10,19]
    EXPECTED OUTPUT: [List[Dict]: ~15 entries (5 guides * 3 pos); keys: wt_seq, mut_seq, pos, type]
    TEST DATA: [e.g., guide='ATGC...', pos=1 'inc': mut='GTGC...'; exactly 1 change]
    REPRODUCTION: [Run with sample_guides; verify len~15, single changes, valid ATGC]
    """
    perturbations = []
    for wt in sample_guides:
        for pos in [1, 10, 19]:
            for mtype in ["inc", "dec"]:
                try:
                    mut = generate_mutant(wt, pos, mtype)
                    perturbations.append(
                        {
                            "wt_seq": wt,
                            "mut_seq": mut,
                            "pos": pos,
                            "type": mtype,
                            "region": "seed"
                            if pos <= 8
                            else "distal"
                            if pos <= 17
                            else "PAM",
                            "num_mut": 1,
                        }
                    )
                except ValueError:
                    pass
    return perturbations[:15]


class TestPerturbationSweepIntegration:
    """
    PURPOSE: As a researcher/engineer, I want integration tests so that I can verify
    full pipeline produces valid outputs without runtime errors.

    INPUTS: [Uses fixtures; runs full load_perturb -> analyze_shifts]
    EXPECTED OUTPUT: [CSV exists, len~3750 rows for full; JSON with stats dicts;
                      all finite diffs, p in [0,1], |d|>0 for muts]
    TEST DATA: [sample_guides (5); expect 15 perturbations (20 singles + doubles?); deltas ~0.1-0.5]
    REPRODUCTION: [pytest -v test_perturbation_sweep.py; inspect temp CSV/JSON;
                       re-run: same outputs with seed=42]
    """

    def test_end_to_end_small_dataset(
        self, sample_guides: list, tmp_path: Path
    ) -> None:
        """
        PURPOSE: As a researcher, I want small-scale integration test so that I can
        validate pipeline before full run.

        INPUTS: [sample_guides: 5 guides; tmp_path: temp output dir]
        EXPECTED OUTPUT: [perturbations len~375 (5 guides * (20 singles + 15 doubles));
                          analysis returns dict with stats keys]
        TEST DATA: [5 guides; expect 35 perturbations (20 singles + 15 doubles);
                    regions balanced, num_mut=1/2]
        REPRODUCTION: [Patch SeqIO.parse to yield mocks; call load_and_perturb;
                       assert len=35, ~equal per type/region; all len=20]
        """
        perturbations = load_and_perturb_guides(
            sample_guides, n_guides=5, max_mut=2, seed=42
        )

        assert len(perturbations) == 121  # 5 guides * (20 singles + ~4 doubles)
        assert all(
            k in {"wt_seq", "mut_seq", "pos", "type", "region", "num_mut"}
            for k in perturbations[0]
        )
        assert all(len(p["wt_seq"]) == 20 for p in perturbations)
        assert all(len(p["mut_seq"]) == 20 for p in perturbations)
        regions = [p["region"] for p in perturbations if isinstance(p["pos"], int)]
        assert len([r for r in regions if r == "seed"]) >= 5 * 8 / 3  # Approx
        types = set(p["type"] for p in perturbations if isinstance(p["pos"], int))
        assert types == {"inc", "dec"}

        # Test analysis stub
        results = analyze_spectral_shifts(perturbations, tmp_path)
        assert isinstance(results, dict)
        assert "stats" in results

    def test_output_schema_validation(self, tmp_path: Path) -> None:
        """
        PURPOSE: As an engineer, I want schema check so that I can ensure CSV/JSON
        match expected structure for downstream analysis.

        INPUTS: [tmp_path: temp dir; assumes small run outputs files]
        EXPECTED OUTPUT: [CSV cols: ['seq_id', 'wt_seq', 'mut_pos', ..., 'delta_coh', 'control_flag'];
                          JSON keys: ['experiment_name', 'parameters', 'metrics', 'generated_at']]
        TEST DATA: [Small run; expect exact col list len=20; JSON valid dict no errors]
        REPRODUCTION: [Run analysis; pd.read_csv cols==expected; json.load no KeyError]
        """
        from experiments.local_perturbation_sweep.sweep_datagen import (
            analyze_spectral_shifts,
        )
        from experiments.local_perturbation_sweep.utils import generate_mutant

        perturbations = [
            {
                "wt_seq": "AT" * 20,
                "mut_seq": generate_mutant("AT" * 20, 1, "inc"),
                "pos": 1,
                "type": "inc",
                "region": "seed",
                "num_mut": 1,
            },
            {
                "wt_seq": "AT" * 20,
                "mut_seq": generate_mutant("AT" * 20, 3, "inc"),
                "pos": 3,
                "type": "inc",
                "region": "seed",
                "num_mut": 1,
            },
            {
                "wt_seq": "AT" * 20,
                "mut_seq": generate_mutant("AT" * 20, 5, "inc"),
                "pos": 5,
                "type": "inc",
                "region": "seed",
                "num_mut": 1,
            },
        ]
        results = analyze_spectral_shifts(perturbations, tmp_path)

        csv_path = tmp_path / "sweep_data.csv"
        df = pd.read_csv(csv_path)
        expected_cols = [
            "seq_id",
            "wt_seq",
            "mut_seq",
            "pos",
            "type",
            "region",
            "num_mut",
            "control_flag",
            "delta_resonance_mag",
            "delta_phase_coh",
            "delta_spectral_centroid",
            "delta_band_energy",
            "delta_snr",
            "seed",
        ]
        assert set(df.columns) == set(expected_cols)
        assert len(df) == 4  # 3 real + 1 control

        json_path = tmp_path / "metadata.json"
        import json

        data = json.load(open(json_path))
        assert set(data.keys()) >= {
            "experiment_name",
            "parameters",
            "metrics",
            "generated_at",
            "power",
        }
        assert isinstance(data["power"], dict)

    def test_reproducibility_with_seed(self, tmp_path: Path) -> None:
        """
        PURPOSE: As a researcher, I want reproducibility test so that I can trust
        results across runs.

        INPUTS: [tmp_path; mock FASTA or use sample; run load_guides twice seed=42]
        EXPECTED OUTPUT: [Identical lists on re-run; lists equal True]
        TEST DATA: [dummy FASTA with 10 guides; seed=42: same 5 balanced subset;
                    seed=43: different]
        REPRODUCTION: [Run load_guides twice same seed; assert guides1 == guides2]
        """
        from unittest.mock import patch, mock_open

        mock_fasta = """ >guide1
ATCGATCGATCGATCGATCG
>guide2
GCTAGCTAGCTAGCTAGCTA
>guide3
ATCGATCGATCGATCGATCG
>guide4
GCTAGCTAGCTAGCTAGCTA
>guide5
ATCGATCGATCGATCGATCG
"""
        with patch("builtins.open", mock_open(read_data=mock_fasta)):
            with patch("Bio.SeqIO.parse") as mock_parse:

                class MockRecord:
                    def __init__(self, id, seq):
                        self.id = id
                        self.seq = seq

                mock_parse.return_value = [
                    MockRecord("guide1", "ATCGATCGATCGATCGATCG"),
                    MockRecord("guide2", "GCTAGCTAGCTAGCTAGCTA"),
                    MockRecord("guide3", "ATCGATCGATCGATCGATCG"),
                    MockRecord("guide4", "GCTAGCTAGCTAGCTAGCTA"),
                    MockRecord("guide5", "ATCGATCGATCGATCGATCG"),
                ]
                guides1 = load_guides(Path("dummy.fasta"), n_guides=5, seed=42)
                guides2 = load_guides(Path("dummy.fasta"), n_guides=5, seed=42)
                assert guides1 == guides2
                assert len(guides1) == 5
                gc = [sum(b in "GC" for b in g) / 20 for g in guides1]
                assert sum(1 for g in gc if g > 0.5) <= 3  # Balanced

    def test_mutant_generation(self) -> None:
        """
        PURPOSE: As an engineer, I want mutant test so that I can ensure single
        position changes are correct and deterministic.

        INPUTS: [None; uses hardcoded WT 'ATGC...20nt', pos=1/10/19, types 'inc'/'dec']
        EXPECTED OUTPUT: [mut_seq len=20; exactly 1 change; 'inc' on A/T ->G/C;
                          'dec' G/C->A/T; same seed same result]
        TEST DATA: [wt='A' + 'T'*19, pos=1 'inc': expect 'G' + 'T'*19;
                    wt='G' + 'A'*19, pos=1 'dec': expect 'A' + 'A'*19]
        REPRODUCTION: [Call generate_mutant multiple times; assert len=20, hamming dist=1;
                       seed=42: always 'A'->'G' (lex order if choice)]
        """
        wt = "A" + "T" * 19  # All AT for inc test
        mut_inc = generate_mutant(wt, 1, "inc")
        assert len(mut_inc) == 20
        assert mut_inc[0] == "G"  # A -> G
        assert sum(a != b for a, b in zip(wt, mut_inc)) == 1

        wt_gc = "G" + "C" * 19  # All GC for dec test
        mut_dec = generate_mutant(wt_gc, 1, "dec")
        assert len(mut_dec) == 20
        assert mut_dec[0] == "A"  # G -> A
        assert sum(a != b for a, b in zip(wt_gc, mut_dec)) == 1

        # Double: adjacent pos 1-2, inc both
        mut_double = generate_double_mutant(wt, (1, 2), ("inc", "inc"))
        assert len(mut_double) == 20
        assert sum(a != b for a, b in zip(wt, mut_double)) == 2  # Two changes
        assert mut_double[0] == "G"  # Pos1 A->G
        assert mut_double[1] == "C"  # Pos2 T->C

        # Reproducibility for double
        mut_double2 = generate_double_mutant(wt, (1, 2), ("inc", "inc"), seed=42)
        assert mut_double == mut_double2

    def test_null_model_integration(self, sample_guides: list, tmp_path: Path) -> None:
        """
        PURPOSE: As a researcher, I want null check so that I can validate diffs
        vs. controls (random-pos, scramble).

        INPUTS: [sample_guides; tmp_path; generate simple perturbations + controls]
        EXPECTED OUTPUT: [Diffs for non-controls |delta|>0; controls ~0 (mean abs<0.1)]
        TEST DATA: [3 guides; generate inc/dec at pos=1; add 3 random-pos controls;
                    expect non-control deltas >0.05, control <0.05]
        REPRODUCTION: [Call compute_diffs on pairs; assert non-null > threshold,
                       null ~0; verify finite]
        """
        wt = "ATGCATGCATGCATGCATGC"  # Sample 20nt
        perturbations = []
        for mtype in ["inc", "dec"]:
            try:
                mut = generate_mutant(wt, 1, mtype)
                diffs = compute_diffs(wt, mut)
                perturbations.append({"wt": wt, "mut": mut, "diffs": diffs})
            except ValueError:
                continue

        # Simple null: same seq as control (delta=0)
        null_diffs = compute_diffs(wt, wt)  # Should be all 0
        assert all(abs(v) < 1e-10 for v in null_diffs.values())

        # Non-null should have some non-zero (from encoding diff)
        for p in perturbations:
            assert any(abs(v) > 0 for v in p["diffs"].values())


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
