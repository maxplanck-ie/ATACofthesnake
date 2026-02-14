import pytest
import pandas as pd
import yaml

@pytest.fixture
def sspath(tmp_path):
    '''
    Create two temporary samplesheets for testing, and two comparison files.
    '''

    sspath = tmp_path / "samplesheets"
    sspath.mkdir()

    bamdir_correct = sspath / "bamdir_correct"
    bamdir_incorrect = sspath / "bamdir_incorrect"

    bamdir_correct.mkdir()
    bamdir_incorrect.mkdir()

    # Create dummy bam files for correct samplesheet
    for sample in ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6', 'sample7', 'sample8']:
        (bamdir_correct / f"{sample}.bam").touch()
    for sample in ['sample1', 'sample2']:
        (bamdir_incorrect / f"{sample}.bam").touch()


    correct = pd.DataFrame({
        'sample': ['sample1', 'sample2', 'sample3', 'sample4', 'sample5' ,'sample6', 'sample7', 'sample8'],
        'factor1': ['A', 'A', 'A', 'A', 'B', 'B', 'B', 'B'],
        'factor2': ['X', 'Y', 'X', 'Y', 'X', 'Y', 'X', 'Y']
    })

    noreps = pd.DataFrame({
        'sample': ['sample1', 'sample2', 'sample3', 'sample4', 'sample5' ,'sample6', 'sample7', 'sample8'],
        'factor1': ['A', 'A', 'A', 'A', 'A', 'A', 'A', 'B'],
        'factor2': ['X', 'Y', 'X', 'Y', 'X', 'Y', 'X', 'Y']
    })

    correctcomp = {
        'comp1': {
            'group1': {'factor1': 'A'},
            'group2': {'factor1': 'B'},
        }
    }
    incorrectcomp = {
        'comp2': {
            'group1': {'factor3': 'A'},
            'group2': {'factor3': 'B'},
        }
    }

    incorrect = pd.DataFrame({
        'factor1': ['A', 'A', 'B', 'B'],
        'factor2': ['X', 'Y', 'X', 'Y']
    })

    correct.to_csv(sspath / "correct.tsv", sep='\t', index=False)
    noreps.to_csv(sspath / "noreps.tsv", sep='\t', index=False)
    incorrect.to_csv(sspath / "incorrect.tsv", sep='\t', index=False)

    with open(sspath / "correctcomp.yaml", 'w') as f:
        yaml.dump(correctcomp, f, default_flow_style=False, sort_keys=False)
    with open(sspath / "incorrectcomp.yaml", 'w') as f:
        yaml.dump(incorrectcomp, f, default_flow_style=False, sort_keys=False)

    return sspath

@pytest.fixture
def fnapath(tmp_path):
    '''
    Create a temporary fasta file for testing.
    '''
    fnapath = tmp_path / "fasta"
    fnapath.mkdir()

    with open(fnapath / "test.fa", 'w') as f:
        f.write(">chr1\n")
        f.write("ATCGATCGATCG\n")
        f.write(">chr2\n")
        f.write("ATCGATCGATCG\n")
        f.write(">chrM\n")
        f.write("ATCGATCGATCG\n")

    with open(fnapath / "rar_bad.bed", 'w') as f:
        f.write("# This is a comment\n")
        f.write("chr1\t100\t200\n")
        f.write("chr2\t300\t400\n")
    
        with open(fnapath / "rar_good.bed", 'w') as f:
            f.write("# This is a comment\n")
            f.write("chr1\t100\t200\n")
            f.write("chr2\t300\t400\n")
            f.write("chrM\t500\t600\n")

    return fnapath