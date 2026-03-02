from pathlib import Path

import pytest

from aos.preflight import Preflight


class TestPreflight:
    def test_samplesheets(self, sspath):
        # No replicates
        with pytest.raises(AssertionError):
            Preflight.validate_samplesheet(
                sspath / "noreps.tsv", sspath / "correctcomp.yaml"
            )
        # Bam files present
        Preflight.validate_samplesheet(
            sspath / "correct.tsv", sspath / "bamdir_correct"
        )
        # Bam files missing
        with pytest.raises(AssertionError):
            Preflight.validate_samplesheet(
                sspath / "correct.tsv", sspath / "bamdir_incorrect"
            )

    def test_optionalpath(self):
        assert Preflight.optional_paths(None) is None
        assert Preflight.optional_paths("test") == Path("test")

    def test_bamdir_validation(self, bampath):
        # Correct bam directory
        Preflight.validate_samples(bampath / "goodbam")
        # Incorrect bam directory
        with pytest.raises(AssertionError):
            Preflight.validate_samples(bampath / "badbam")
        # Incorrect cram directory
        with pytest.raises(AssertionError):
            Preflight.validate_samples(bampath / "badcram")

    def test_mitostring_validation(self, fnapath):
        # Correct mitostring
        Preflight.validate_mitostring(
            fnapath / "test.fa", "chrM", fnapath / "rar_good.bed"
        )
        # Incorrect mitostring
        with pytest.raises(AssertionError):
            Preflight.validate_mitostring(
                fnapath / "test.fa", "MT", fnapath / "rar_good.bed"
            )
        # Mitostring in rar file
        Preflight.validate_mitostring(
            fnapath / "test.fa", "chrM", fnapath / "rar_good.bed"
        )
        # Mitostring not in rar file
        with pytest.raises(AssertionError):
            Preflight.validate_mitostring(
                fnapath / "test.fa", "chrM", fnapath / "rar_bad.bed"
            )
