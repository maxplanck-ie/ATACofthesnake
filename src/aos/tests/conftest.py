from copy import deepcopy

import pandas as pd
import pytest
import yaml


@pytest.fixture
def sspath(tmp_path):
    """
    Create two temporary samplesheets for testing, and two comparison files.
    """

    sspath = tmp_path / "samplesheets"
    sspath.mkdir()

    bamdir_correct = sspath / "bamdir_correct"
    bamdir_incorrect = sspath / "bamdir_incorrect"

    bamdir_correct.mkdir()
    bamdir_incorrect.mkdir()

    # Create dummy bam files for correct samplesheet
    for sample in [
        "sample1",
        "sample2",
        "sample3",
        "sample4",
        "sample5",
        "sample6",
        "sample7",
        "sample8",
    ]:
        (bamdir_correct / f"{sample}.bam").touch()
    for sample in ["sample1", "sample2"]:
        (bamdir_incorrect / f"{sample}.bam").touch()

    correct = pd.DataFrame(
        {
            "sample": [
                "sample1",
                "sample2",
                "sample3",
                "sample4",
                "sample5",
                "sample6",
                "sample7",
                "sample8",
            ],
            "factor1": ["A", "A", "A", "A", "B", "B", "B", "B"],
            "factor2": ["X", "Y", "X", "Y", "X", "Y", "X", "Y"],
        }
    )

    noreps = pd.DataFrame(
        {
            "sample": [
                "sample1",
                "sample2",
                "sample3",
                "sample4",
                "sample5",
                "sample6",
                "sample7",
                "sample8",
            ],
            "factor1": ["A", "A", "A", "A", "A", "A", "A", "B"],
            "factor2": ["X", "Y", "X", "Y", "X", "Y", "X", "Y"],
        }
    )

    correctcomp = {
        "comp1": {
            "group1": {"factor1": "A"},
            "group2": {"factor1": "B"},
        }
    }
    incorrectcomp = {
        "comp2": {
            "group1": {"factor3": "A"},
            "group2": {"factor3": "B"},
        }
    }

    incorrect = pd.DataFrame(
        {"factor1": ["A", "A", "B", "B"], "factor2": ["X", "Y", "X", "Y"]}
    )

    correct.to_csv(sspath / "correct.tsv", sep="\t", index=False)
    noreps.to_csv(sspath / "noreps.tsv", sep="\t", index=False)
    incorrect.to_csv(sspath / "incorrect.tsv", sep="\t", index=False)

    with open(sspath / "correctcomp.yaml", "w") as f:
        yaml.dump(correctcomp, f, default_flow_style=False, sort_keys=False)
    with open(sspath / "incorrectcomp.yaml", "w") as f:
        yaml.dump(incorrectcomp, f, default_flow_style=False, sort_keys=False)

    return sspath


@pytest.fixture
def fnapath(tmp_path):
    """
    Create a temporary fasta file for testing.
    """
    fnapath = tmp_path / "fasta"
    fnapath.mkdir()

    with open(fnapath / "test.fa", "w") as f:
        f.write(">chr1\n")
        f.write("ATCGATCGATCG\n")
        f.write(">chr2\n")
        f.write("ATCGATCGATCG\n")
        f.write(">chrM\n")
        f.write("ATCGATCGATCG\n")

    with open(fnapath / "rar_bad.bed", "w") as f:
        f.write("# This is a comment\n")
        f.write("chr1\t100\t200\n")
        f.write("chr2\t300\t400\n")

        with open(fnapath / "rar_good.bed", "w") as f:
            f.write("# This is a comment\n")
            f.write("chr1\t100\t200\n")
            f.write("chr2\t300\t400\n")
            f.write("chrM\t500\t600\n")

    return fnapath


@pytest.fixture
def bampath(tmp_path):
    """
    Create a directory with bamfiles for testing
    """
    bampath = tmp_path / "bams"
    bampath.mkdir()

    goodbam = bampath / "goodbam"
    badbam = bampath / "badbam"
    badcram = bampath / "badcram"
    goodbam.mkdir()
    badbam.mkdir()
    badcram.mkdir()

    (goodbam / "file1.bam").touch()
    (goodbam / "file2.bam").touch()
    (goodbam / "file3.cram").touch()

    (badbam / "file-1.bam").touch()
    (badcram / "file-1.cram").touch()

    return bampath


@pytest.fixture
def comppath(tmp_path):
    """
    Create a samplesheet and comparison file for comparison validation testing
    """

    def wy(d, f):
        with open(f, "w") as of:
            yaml.dump(d, of)

    with open(tmp_path / "ss.tsv", "w") as f:
        f.write("sample\ttime\ttreatment\n")
        f.write("sample1\t0\tcontrol\n")
        f.write("sample2\t0\tcontrol\n")
        f.write("sample3\t0\ttreatment\n")
        f.write("sample4\t0\ttreatment\n")
        f.write("sample5\t1\tcontrol\n")
        f.write("sample6\t1\tcontrol\n")
        f.write("sample7\t1\ttreatment\n")
        f.write("sample8\t1\ttreatment\n")
        f.write("sample9\t2\tcontrol\n")
        f.write("sample10\t2\tcontrol\n")
        f.write("sample11\t2\ttreatment\n")
        f.write("sample12\t2\ttreatment\n")

    # correct comparison
    comp = {
        "test_comp": {
            "type": "twogroup",
            "design": "~time",
            "t0_1": {"time": [0, 1]},
            "t2": {"time": 2},
        }
    }
    wy(comp, tmp_path / "comp.yaml")
    # notype
    d = deepcopy(comp)
    d["test_comp"].pop("type")
    wy(d, tmp_path / "comp_notype.yaml")
    # wrongtype
    d = deepcopy(comp)
    d["test_comp"]["type"] = "wrongtype"
    wy(d, tmp_path / "comp_wrongtype.yaml")
    # notilde
    d = deepcopy(comp)
    d["test_comp"]["design"] = "time"
    wy(d, tmp_path / "comp_notilde.yaml")
    # var not in ss
    d = deepcopy(comp)
    d["test_comp"]["design"] = "~time + factortype"
    wy(d, tmp_path / "comp_varnotins.yaml")
    # principle of marginality
    d = deepcopy(comp)
    d["test_comp"]["design"] = "~time + time:treatment"
    wy(d, tmp_path / "comp_margprin.yaml")

    # tg - not2
    d = deepcopy(comp)
    d["test_comp"].pop("t2")
    wy(d, tmp_path / "comp_tg_not2.yaml")
    # tg - wrongfactor
    d = deepcopy(comp)
    d["test_comp"]["t2"] = {"blime": [2]}
    wy(d, tmp_path / "comp_tg_wrongfactor.yaml")
    # tg - wronglevel
    d = deepcopy(comp)
    d["test_comp"]["t2"] = {"time": 3}
    wy(d, tmp_path / "comp_tg_wronglevel.yaml")
    # tg - wronglevel
    d = deepcopy(comp)
    d["test_comp"]["t0_1"] = [0, 3]
    wy(d, tmp_path / "comp_tg_wronglevel2.yaml")

    # lrt
    comp = {
        "test_comp": {"type": "lrt", "design": "~time*treatment", "reduced": "~time"}
    }
    wy(comp, tmp_path / "comp_lrt.yaml")
    # lrt - noreduced
    d = deepcopy(comp)
    d["test_comp"].pop("reduced")
    wy(d, tmp_path / "comp_lrt_noreduced.yaml")
    # lrt - notilde
    d = deepcopy(comp)
    d["test_comp"]["reduced"] = "time"
    wy(d, tmp_path / "comp_lrt_notilde.yaml")
    # lrt - wrong reduced factor
    d = deepcopy(comp)
    d["test_comp"]["reduced"] = "~blime"
    wy(d, tmp_path / "comp_lrt_wrongreduced.yaml")
    # lrt - reduced not reduced
    d = deepcopy(comp)
    d["test_comp"]["reduced"] = "~time*treatment"
    wy(d, tmp_path / "comp_lrt_reducednotreduced.yaml")
    # lrt - notnested
    d = deepcopy(comp)
    d["test_comp"]["design"] = "~time+time:treatment"
    d["test_comp"]["reduced"] = "~treatment+time:treatment"
    wy(d, tmp_path / "comp_lrt_notnested.yaml")

    # gp
    comp = {
        "test_comp": {"type": "timecourse", "time": "time", "time_type": "continuous"}
    }
    wy(comp, tmp_path / "comp_gp.yaml")

    # gp - notime
    d = deepcopy(comp)
    d["test_comp"].pop("time")
    wy(d, tmp_path / "comp_gp_notime.yaml")
    # gp - notime_type
    d = deepcopy(comp)
    d["test_comp"].pop("time_type")
    wy(d, tmp_path / "comp_gp_notime_type.yaml")
    # gp - wrong time_type
    d = deepcopy(comp)
    d["test_comp"]["time_type"] = "wrong"
    wy(d, tmp_path / "comp_gp_wrongtime_type.yaml")
    # gp - wrong time
    d = deepcopy(comp)
    d["test_comp"]["time"] = "wrong"
    wy(d, tmp_path / "comp_gp_wrongtime.yaml")

    # gp - ordinal
    comp = {
        "test_comp": {
            "type": "timecourse",
            "time": "time",
            "time_type": "ordinal",
            "order": [0, 1, 2],
        }
    }
    wy(comp, tmp_path / "comp_gp_ordinal.yaml")
    # gp - noorder
    d = deepcopy(comp)
    d["test_comp"].pop("order")
    wy(d, tmp_path / "comp_gp_noorder.yaml")
    # gp - order wrong level
    d = deepcopy(comp)
    d["test_comp"]["order"] = [0, 1, 5]
    wy(d, tmp_path / "comp_gp_ord_wronglevel.yaml")

    return tmp_path
