import pytest
import pandas as pd
from saga.plot import get_bin_kinematics_title, get_lims_coords, get_bin_centers


@pytest.fixture(name="df")
def df_fixture():
    return pd.DataFrame(
        {
            "x": [0.1, 0.3],
            "Q2": [1.0, 3.0],
            "x_err": [0.05, 0.05],
            "Q2_err": [0.05, 0.05],
        }
    )


@pytest.fixture(name="cols")
def cols_fixture():
    return ["x", "Q2"]


@pytest.fixture(name="col_titles")
def col_titles_fixture():
    return {"x": "x", "Q2": "Q^{2} (GeV^{2})"}


@pytest.fixture(name="node")
def node_fixture():
    return {
        "x": [0.0, 0.5, 1.0],
        "Q2": [0.0, 5.0, 11.0],
    }


@pytest.fixture(name="node_nested")
def node_nested_fixture():
    return {
        "nested": [
            {
                "Q2": {
                    "nbins": 2,
                    "lims": [0.0, 5.0, 11.0],
                    "nested": [
                        {
                            "x": {
                                "nbins": 2,
                                "lims": [0.0, 0.5, 1.0],
                            }
                        },
                        {
                            "x": {
                                "nbins": 2,
                                "lims": [0.0, 0.5, 1.0],
                            }
                        },
                    ],
                }
            },
        ]
    }


@pytest.fixture(name="outer_xlims")
def outer_xlims_fixture():
    return [0.0, 1.0]


@pytest.fixture(name="outer_ylims")
def outer_ylims_fixture():
    return [0.0, 11.0]


@pytest.fixture(name="cuts")
def cuts_fixture():
    return {
        0: "(x>=0.0 && x<0.5) && (Q2>=0.0 && Q2<5.0)",
        1: "(x>=0.0 && x<0.5) && (Q2>=5.0 && Q2<11.0)",
        2: "(x>=0.5 && x<1.0) && (Q2>=0.0 && Q2<5.0)",
        3: "(x>=0.5 && x<1.0) && (Q2>=5.0 && Q2<11.0)",
    }


def test_get_bin_kinematics_title(df, cols, col_titles):

    # Test first and second bin titles
    assert (
        get_bin_kinematics_title(0, df, cols=cols, col_titles=col_titles)
        == "$<x> = 0.10\\pm0.05$ , $<Q^{2} (GeV^{2})> = 1.00\\pm0.05$"
    )
    assert (
        get_bin_kinematics_title(1, df, cols=cols, col_titles=col_titles)
        == "$<x> = 0.30\\pm0.05$ , $<Q^{2} (GeV^{2})> = 3.00\\pm0.05$"
    )


def test_get_lims_coords(node, node_nested, outer_xlims, outer_ylims):

    # Get grid limits coordinates
    lims_coords_grid = get_lims_coords(
        node,
        outer_xlims,
        outer_ylims,
        var_keys=["x", "Q2"],
        nested_key="nested",
        lims_key="lims",
        swap_axes=False,
    )
    assert lims_coords_grid == [[[0.5, 0.5], [0.0, 11.0]], [[0.0, 1.0], [5.0, 5.0]]]

    # Get nested limits coordinates
    lims_coords_nested = get_lims_coords(
        node_nested,
        outer_xlims,
        outer_ylims,
        var_keys=None,
        nested_key="nested",
        lims_key="lims",
        swap_axes=False,
    )
    assert lims_coords_nested == [
        [[0.5, 0.5], [0.0, 5.0]],
        [[0.5, 0.5], [5.0, 11.0]],
        [[0.0, 1.0], [5.0, 5.0]],
    ]

    # Test errors
    with pytest.raises(ValueError):
        get_lims_coords({}, outer_xlims, outer_ylims)


def test_get_bin_centers(cuts):

    # Test centers and widths
    bin_centers, bin_widths = get_bin_centers(cuts)
    assert bin_centers == {
        0: [0.25, 2.5],
        1: [0.25, 8.0],
        2: [0.75, 2.5],
        3: [0.75, 8.0],
    }
    assert bin_widths == {
        0: [0.5, 5.0],
        1: [0.5, 6.0],
        2: [0.5, 5.0],
        3: [0.5, 6.0],
    }
