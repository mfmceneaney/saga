import pytest
from saga.aggregate import (
    get_config_str,
    get_config_out_path,
    get_out_file_name,
    get_config_list,
    get_out_dirs_list,
)


@pytest.fixture(name="config_simple")
def config_simple_fixture():
    return {
        "string": "abcde",
        "float": 0.0,
        "int": 0,
        "str_list": ["a", "b", "c", "d", "e"],
        "float_list": [0.0, 1.0, 2.0, 3.0, 4.0],
        "int_list": [0, 1, 2, 3, 4],
    }


@pytest.fixture(name="config_complex")
def config_complex_fixture():
    return {
        "binschemes": {"binscheme1": {"binvar1": [0.0, 1.0], "binvar2": [0.0, 2.0]}}
    }


@pytest.fixture(name="sep")
def sep_fixture():
    return "_"


@pytest.fixture(name="aliases")
def aliases_fixture():
    return {
        "binschemes": {
            "{'binscheme1': {'binvar1': [0.0, 1.0], 'binvar2': [0.0, 2.0]}}": "binscheme1_value"
        }
    }


@pytest.fixture(name="base_dir")
def base_dir_fixture():
    return "/path/to/base_dir"


@pytest.fixture(name="aggregate_keys")
def aggregate_keys_fixture():
    return ["string"]


@pytest.fixture(name="result_name")
def result_name_fixture():
    return "a0"


@pytest.fixture(name="ext")
def ext_fixture():
    return ".csv"


@pytest.fixture(name="configs")
def configs_fixture():
    return {"string": ["abcde", "fghij"], "float": [0.0, 1.0]}


@pytest.fixture(name="out_dirs_list")
def out_dirs_list_fixture():
    return [
        [
            "/path/to/base_dir/float_0.0__string_abcde",
            "/path/to/base_dir/float_0.0__string_fghij",
        ],
        [
            "/path/to/base_dir/float_1.0__string_abcde",
            "/path/to/base_dir/float_1.0__string_fghij",
        ],
    ]


def test_get_config_str(config_simple, config_complex, sep, aliases):

    # Test str, float and int cases
    config_str = get_config_str(config_simple, sep=sep, aliases={})
    assert (
        config_str
        == "float_0.0__float_list_0.0_1.0_2.0_3.0_4.0__int_0__int_list_0_1_2_3_4__str_list_a_b_c_d_e__string_abcde"
    )

    # Test alias case
    config_str = get_config_str(config_complex, sep=sep, aliases=aliases)
    assert config_str == "binscheme1_value"


def test_get_config_out_path(
    base_dir,
    aggregate_keys,
    result_name,
    config_simple,
    config_complex,
    sep,
    aliases,
    ext,
):

    # Test simple config without aggregate keys
    config_out_path = get_config_out_path(
        base_dir, [], result_name, config_simple, sep=sep, aliases={}, ext=ext
    )
    assert (
        config_out_path
        == "/path/to/base_dir/aggregate______float_0.0__float_list_0.0_1.0_2.0_3.0_4.0__int_0__int_list_0_1_2_3_4__str_list_a_b_c_d_e__string_abcde___a0.csv"
    )

    # Test simple config with aggregate keys
    config_out_path = get_config_out_path(
        base_dir,
        aggregate_keys,
        result_name,
        config_simple,
        sep=sep,
        aliases={},
        ext=ext,
    )
    assert (
        config_out_path
        == "/path/to/base_dir/aggregate___string___float_0.0__float_list_0.0_1.0_2.0_3.0_4.0__int_0__int_list_0_1_2_3_4__str_list_a_b_c_d_e__string_abcde___a0.csv"
    )

    # Test complex config without aggregate keys
    config_out_path = get_config_out_path(
        base_dir, [], result_name, config_complex, sep=sep, aliases=aliases, ext=ext
    )
    assert (
        config_out_path == "/path/to/base_dir/aggregate______binscheme1_value___a0.csv"
    )

    # Test complex config with aggregate keys
    config_out_path = get_config_out_path(
        base_dir,
        aggregate_keys,
        result_name,
        config_complex,
        sep=sep,
        aliases=aliases,
        ext=ext,
    )
    assert (
        config_out_path
        == "/path/to/base_dir/aggregate___string___binscheme1_value___a0.csv"
    )


def test_get_out_file_name(base_dir, ext):

    # Test simple name formation
    assert get_out_file_name(None, "", "x", ext=ext) == "x.csv"
    assert get_out_file_name(None, "base_", "x", ext=ext) == "base_x.csv"
    assert (
        get_out_file_name(base_dir, "base_", "x", ext=ext)
        == "/path/to/base_dir/base_x.csv"
    )


def test_get_config_list(configs, aggregate_keys):

    # Test without aggregate keys
    assert get_config_list(configs) == [
        {"string": "abcde", "float": 0.0},
        {"string": "abcde", "float": 1.0},
        {"string": "fghij", "float": 0.0},
        {"string": "fghij", "float": 1.0},
    ]

    # And with some aggregate keys
    assert get_config_list(configs, aggregate_keys=aggregate_keys) == [
        {"float": 0.0},
        {"float": 1.0},
    ]

    # And with all aggregate keys
    assert get_config_list(configs, aggregate_keys=list(configs.keys())) == []


def test_get_out_dirs_list(configs, base_dir, aggregate_keys, aliases, out_dirs_list):

    # Test simple configs
    assert (
        get_out_dirs_list(
            configs, base_dir, aggregate_keys=aggregate_keys, aliases=aliases
        )
        == out_dirs_list
    )
