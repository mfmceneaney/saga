import pytest
import os
import yaml
from saga.orchestrate import create_jobs


@pytest.fixture(name="configs")
def configs_fixture():
    return {"string": ["abcde", "fghij"], "float": [0.0, 1.0]}


@pytest.fixture(name="yaml_args")
def yaml_args_fixture():
    return {"float": 0.0}


@pytest.fixture(name="out_dirs_list")
def out_dirs_list_fixture():
    return [
        "float_0.0__string_abcde",
        "float_0.0__string_fghij",
        "float_1.0__string_abcde",
        "float_1.0__string_fghij",
    ]


@pytest.fixture(name="out_configs")
def out_configs_fixture():
    return [
        {"float": 0.0, "string": "abcde"},
        {"float": 0.0, "string": "fghij"},
        {"float": 1.0, "string": "abcde"},
        {"float": 1.0, "string": "fghij"},
    ]


def test_create_jobs(
    tmp_path,
    configs,
    yaml_args,
    out_dirs_list,
    out_configs,
):

    # Set some local variables
    submit_file_name = "submit.sh"
    submit_file_content = "#!/bin/bash\n"
    yaml_file_name = "args.yaml"

    # Write submit script
    submit_path = os.path.join(tmp_path, submit_file_name)
    with open(submit_path, "w", encoding="utf-8") as f:
        f.write(submit_file_content)

    # Write yaml file
    yaml_path = os.path.join(tmp_path, yaml_file_name)
    with open(yaml_path, "w", encoding="utf-8") as f:
        yaml.safe_dump(yaml_args, f, default_flow_style=False)

    # Create simple config directories
    create_jobs(
        configs, tmp_path, submit_path=submit_path, yaml_path=yaml_path, aliases=None
    )

    # Loop config directories
    for idx, path_i in enumerate(out_dirs_list):

        # Check existence
        assert os.path.isdir(os.path.join(tmp_path, path_i))
        submit_path_i = os.path.join(tmp_path, path_i, submit_file_name)
        assert os.path.isfile(submit_path_i)
        yaml_path_i = os.path.join(tmp_path, path_i, yaml_file_name)
        assert os.path.isfile(yaml_path_i)

        # Get yaml args for config
        yaml_args_i = out_configs[idx]
        yaml_args_i["outdir"] = os.path.join(tmp_path, path_i)

        # Check file contents
        with open(submit_path_i, "r", encoding="utf-8") as f:
            assert f.read() == submit_file_content
        with open(yaml_path_i, "r", encoding="utf-8") as f:
            assert yaml.safe_load(f) == yaml_args_i
