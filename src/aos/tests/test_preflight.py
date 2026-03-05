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

    def test_comparison_validation(self, comppath):
        ss = comppath / "ss.tsv"
        # Correct comp
        Preflight.validate_comparison(comppath / "comp.yaml", ss)
        # no type
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_notype.yaml", ss)
        # wrong type
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_wrongtype.yaml", ss)
        # notilde
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_notilde.yaml", ss)
        # var not in ss
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_varnotins.yaml", ss)
        # Marginality principle
        _ = Preflight.validate_comparison(comppath / "comp_margprin.yaml", ss)
        assert len(_) == 1

    def test_comparison_twogroup(self, comppath):
        ss = comppath / "ss.tsv"
        # Not two groups
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_tg_not2.yaml", ss)
        # wrong factor (not in ss)
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_tg_wrongfactor.yaml", ss)
        # wrong level (not in ss)
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_tg_wronglevel.yaml", ss)
        # wrong level (not in ss) for list
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_tg_wronglevel2.yaml", ss)

    def test_comparison_lrt(self, comppath):
        ss = comppath / "ss.tsv"
        # lrt
        Preflight.validate_comparison(comppath / "comp_lrt.yaml", ss)
        # lrt - noreduced
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_lrt_noreduced.yaml", ss)
        # lrt - notilde
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_lrt_notilde.yaml", ss)
        # lrt - wrong reduced factor
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_lrt_wrongreduced.yaml", ss)
        # lrt - reduced not reduced
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(
                comppath / "comp_lrt_reducednotreduced.yaml", ss
            )
        # lrt - notnested
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_lrt_notnested.yaml", ss)

    def test_comparison_gp(self, comppath):
        ss = comppath / "ss.tsv"
        # gp
        Preflight.validate_comparison(comppath / "comp_gp.yaml", ss)
        # gp - notime
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_gp_notime.yaml", ss)
        # gp - notime_type
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_gp_notime_type.yaml", ss)
        # gp - wrong type
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_gp_wrongtime_type.yaml", ss)
        # gp - wrong time
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_gp_wrongtime.yaml", ss)
        # gp - ordinal
        Preflight.validate_comparison(comppath / "comp_gp_ordinal.yaml", ss)
        # gp - noorder
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_gp_noorder.yaml", ss)
        # gp - order wrong level
        with pytest.raises(AssertionError):
            Preflight.validate_comparison(comppath / "comp_gp_ord_wronglevel.yaml", ss)
