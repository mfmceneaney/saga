import pytest
import yaml
import numpy as np
import pandas as pd
import uproot as ur
from saga.data import load_th1, load_yaml, load_csv, save_txt

@pytest.fixture
def config_simple():
    return {
        'string': 'abcde',
        'float': 0.0,
        'int': 0,
        'str_list': ['a', 'b', 'c', 'd', 'e'],
        'float_list': [0.0, 1.0, 2.0, 3.0, 4.0],
        'int_list': [0, 1, 2, 3, 4]
    }

@pytest.fixture
def config_complex():
    return {
        "binschemes":{
            "binscheme1":{
                "binvar1":[0.0,1.0],
                "binvar2":[0.0,2.0]
            }
        }
    }

@pytest.fixture
def sep():
    return "_"

@pytest.fixture
def aliases():
    return {
        'binschemes': {
            "{'binscheme1': {'binvar1': [0.0, 1.0], 'binvar2': [0.0, 2.0]}}": 'binscheme1_value'
        }
    }

@pytest.fixture
def base_dir():
    return "/path/to/base_dir"

@pytest.fixture
def aggregate_keys():
    return ["string"]

@pytest.fixture
def result_name():
    return "a0"

@pytest.fixture
def ext():
    return ".csv"

@pytest.fixture
def configs():
    return {'string': ['abcde', 'fghij'], 'float': [0.0, 1.0]}

@pytest.fixture
def out_dirs_list():
    return [
        [
            '/path/to/base_dir/float_0.0__string_abcde',
            '/path/to/base_dir/float_0.0__string_fghij',
        ],
        [
            '/path/to/base_dir/float_1.0__string_abcde',
            '/path/to/base_dir/float_1.0__string_fghij',
        ],
    ]

def test_load_th1():

    # Set function parameters
    path = "test_th1.root"
    name = "th1"

    # Define histogram parameters
    n_bins = 2
    xlims = (-1,1)
    bin_edges = np.linspace(*xlims, n_bins + 1)

    # Set some dummy data
    x = [-0.5,]

    # Create the numpy histogram
    hist = np.histogram(x, bins=bin_edges)

    # Write to ROOT file using uproot
    with ur.recreate(path) as f:
        f[name] = hist

    # Test histogram loading
    th1 = load_th1(path,name=name)
    assert np.all(th1[0] == [1., 0.]) and np.all(th1[1] == [-1., 0., 1.])

    # Test default histogram loading
    assert load_th1("someotherhistogram.root",name=name) == []

def test_load_yaml(config_simple):

    # Write config to YAML
    yaml_path = "test.yaml"
    with open(yaml_path, 'w') as yaml_file:
            yaml.safe_dump(config_simple, yaml_file, default_flow_style=False)

    # Load back from file and check
    assert load_yaml(yaml_path) == config_simple

def test_load_csv(config_simple,config_complex,aggregate_keys,aliases):

    # Write simple csv
    path = "test.csv"
    data = np.array([[0.1,0.2],[-0.1,-0.2]],dtype=float)
    data_pd = pd.DataFrame(data)
    save_txt(path,data)
    loaded_pd = load_csv(path)
    assert np.array_equal(loaded_pd.values, data_pd.values)

#TODO: Expand testing above and testing for additional methods below.
