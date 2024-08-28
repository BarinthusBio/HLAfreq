import pandas as pd
import numpy as np
import pytest
import HLAfreq
import HLAfreq.HLAfreq_pymc as HLAhdi


@pytest.fixture(
    params=[
        ["A*01:01", "A*01:02", "A*01:03"],
        ["A*01:01:01", "A*01:02:01", "A*01:03:01"],
        ["A*01:01", "A*01:02", "A*01:03:01"],
        ["A*01:01", "A*01:02", "A*01:03:01", "A*01:03:02"],
    ]
)
def alleles(request):
    return request.param


@pytest.fixture
def aftab(alleles):
    return HLAfreq.simulate_study(alleles, 3, "X")


# Not a meaningful test but slow
# def test_hdi():
#     hdi = HLAhdi.AFhdi(aftab, credible_interval=0.95)
#     assert all(hdi.columns == ["lo", "hi", "allele", "post_mean"])


def test_correct_c_array(aftab):
    aftab = HLAfreq.only_complete(aftab)
    aftab = HLAfreq.decrease_resolution(aftab, 2)
    aftab['c'] = 2 * aftab.allele_freq * aftab.sample_size
    c_pivot = aftab.pivot(columns="allele", index="population", values="c")
    c_array = HLAhdi._make_c_array(aftab)
    pytest.approx(c_array[0]) == c_pivot


def test_correct_c_array_alleles(aftab):
    aftab = HLAfreq.only_complete(aftab)
    aftab = HLAfreq.decrease_resolution(aftab, 2)
    aftab['c'] = 2 * aftab.allele_freq * aftab.sample_size
    c_pivot = aftab.pivot(columns="allele", index="population", values="c")
    c_array = HLAhdi._make_c_array(aftab)
    c_array[1] == c_pivot.columns


def test_complete(aftab):
    aftab = aftab.drop(0)
    with pytest.raises(AssertionError) as e:
        HLAhdi._make_c_array(aftab)
    assert "incomplete_studies" in str(e.value)


def test_unique_in_study(aftab):
    aftab["allele"] = "A*00:00"
    with pytest.raises(AssertionError) as e:
        caf = HLAhdi._make_c_array(aftab)
    assert "same allele appears multiple times" in str(e.value)


def test_resolution_check(aftab):
    aftab.loc[0, "allele"] = aftab.loc[0, "allele"] + ":00:00"
    with pytest.raises(AssertionError) as e:
        caf = HLAhdi._make_c_array(aftab)
    assert "multiple resolutions" in str(e.value)

