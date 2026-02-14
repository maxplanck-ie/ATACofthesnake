import pytest
from aos.preflight import Preflight
from pathlib import Path

class TestPreflight:
    def test_samplesheets(self, sspath):
        # No replicates
        with pytest.raises(AssertionError):
            Preflight.validate_samplesheet(sspath / "noreps.tsv", sspath / "correctcomp.yaml" )
        # Bam files present
        Preflight.validate_samplesheet(sspath / "correct.tsv", sspath / "bamdir_correct")
        # Bam files missing
        with pytest.raises(AssertionError):
            Preflight.validate_samplesheet(sspath / "correct.tsv", sspath / "bamdir_incorrect")
    
    # def test_comparison(self, sspath):
    #     # Both correct
    #     Preflight.validate_comparison(sspath / "correct.tsv", sspath / "correctcomp.yaml" )
    #     # Wrong factor in comparison file
    #     with pytest.raises(AssertionError):
    #         Preflight.validate_comparison(sspath / "correct.tsv", sspath / "incorrectcomp.yaml" )
    #     # No sample
    #     with pytest.raises(AssertionError):
    #         Preflight.validate_comparison(sspath / "incorrect.tsv", sspath / "correctcomp.yaml" )
    
    def test_optionalpath(self):
        assert Preflight.optional_paths(None) is None
        assert Preflight.optional_paths("test") == Path("test")
    
    def test_mitostring_validation(self, fnapath):
        # Correct mitostring
        Preflight.validate_mitostring(fnapath / "test.fa", "chrM", fnapath / "rar_good.bed")
        # Incorrect mitostring
        with pytest.raises(AssertionError):
            Preflight.validate_mitostring(fnapath / "test.fa", "MT", fnapath / "rar_good.bed")
        # Mitostring in rar file
        Preflight.validate_mitostring(fnapath / "test.fa", "chrM", fnapath / "rar_good.bed")
        # Mitostring not in rar file
        with pytest.raises(AssertionError):
            Preflight.validate_mitostring(fnapath / "test.fa", "chrM", fnapath / "rar_bad.bed")
        