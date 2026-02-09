import pytest
from aos.preflight import Preflight
from pathlib import Path

class TestPreflight:
    def test_samplesheets(self, sspath):
        # Both correct
        Preflight.validate_comparison(sspath / "correct.tsv", sspath / "correctcomp.yaml" )
        # Wrong factor in comparison file
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(sspath / "correct.tsv", sspath / "incorrectcomp.yaml" )
        # No sample
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(sspath / "incorrect.tsv", sspath / "correctcomp.yaml" )
    
    def test_optionalpath(self):
        assert Preflight.optional_paths(None) is None
        assert Preflight.optional_paths("test") == Path("test")
    
    def test_mitostring_validation(self, fnapath):
        # Correct mitostring
        Preflight.validate_mitostring(fnapath / "test.fa", "chrM")
        # Incorrect mitostring
        with pytest.raises(AssertionError):
            Preflight.validate_mitostring(fnapath / "test.fa", "MT")