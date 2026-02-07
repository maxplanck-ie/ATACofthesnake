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

    correct = pd.DataFrame({
        'sample': ['sample1', 'sample2', 'sample3', 'sample4'],
        'factor1': ['A', 'A', 'B', 'B'],
        'factor2': ['X', 'Y', 'X', 'Y']
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
    incorrect.to_csv(sspath / "incorrect.tsv", sep='\t', index=False)

    with open(sspath / "correctcomp.yaml", 'w') as f:
        yaml.dump(correctcomp, f, default_flow_style=False, sort_keys=False)
    with open(sspath / "incorrectcomp.yaml", 'w') as f:
        yaml.dump(incorrectcomp, f, default_flow_style=False, sort_keys=False)

    return sspath