import pytest
from AOS.preflight import Preflight
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
        assert Preflight.optional_paths(None) == None
        assert Preflight.optional_paths("test") == Path("test")