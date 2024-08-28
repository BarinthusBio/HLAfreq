import pytest
import HLAfreq
import pandas as pd
import numpy as np
from collections.abc import Iterable


@pytest.fixture(params=["Thailand"])
def country(request):
    return request.param


@pytest.fixture
def base_url(country):
    return HLAfreq.makeURL(country)


def test_url_string(base_url):
    assert isinstance(base_url, str)


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


@pytest.fixture
def expected_freq(aftab):
    aftab = HLAfreq.only_complete(aftab)
    aftab = HLAfreq.decrease_resolution(aftab, 2)
    aftab = HLAfreq.unmeasured_alleles(aftab, "population")
    grouped = aftab.groupby("allele")[["allele_freq", "sample_size"]]
    return grouped.apply(
        lambda x: np.average(x.allele_freq, weights=x.sample_size * 2)
    ).values


def test_wav(aftab, expected_freq):
    aftab = HLAfreq.only_complete(aftab)
    aftab = HLAfreq.decrease_resolution(aftab, 2)
    caf = HLAfreq.combineAF(aftab)
    assert caf.wav.values == pytest.approx(expected_freq)


def test_combineAF(aftab, expected_freq):
    aftab = HLAfreq.only_complete(aftab)
    aftab = HLAfreq.decrease_resolution(aftab, 2)
    caf = HLAfreq.combineAF(aftab)
    assert caf.allele_freq.values == pytest.approx(expected_freq, rel=5e-2)


def test_resolution_check(aftab):
    aftab.loc[0, "allele"] = aftab.loc[0, "allele"] + ":00:00"
    with pytest.raises(AssertionError) as e:
        caf = HLAfreq.combineAF(aftab)
    assert "multiple resolutions" in str(e.value)


def test_low_res(aftab):
    with pytest.raises(AssertionError) as e:
        HLAfreq.decrease_resolution(aftab, 5)
    assert "resolution below" in str(e.value)


def test_single_locus(aftab):
    with pytest.raises(AssertionError) as e:
        aftab.iloc[0, 1] = "Y"
        HLAfreq.combineAF(aftab)
    assert "only 1 loci" in str(e.value)


def test_complete(aftab):
    aftab = aftab.drop(0)
    with pytest.raises(AssertionError) as e:
        HLAfreq.combineAF(aftab)
    assert "incomplete_studies" in str(e.value)


def test_single_sample_size(aftab):
    aftab = HLAfreq.only_complete(aftab)
    aftab = HLAfreq.decrease_resolution(aftab, 2)
    aftab.loc[0, "sample_size"] = 51
    with pytest.raises(AssertionError) as e:
        caf = HLAfreq.combineAF(aftab)
    assert "dataset_sample_size must be 1" in str(e.value)


def test_unique_in_study(aftab):
    aftab["allele"] = "A*00:00"
    with pytest.raises(AssertionError) as e:
        caf = HLAfreq.combineAF(aftab)
    assert "same allele appears multiple times" in str(e.value)


def test_multi_loci_allele_decrease_resolution(aftab):
    aftabB = aftab.copy()
    aftabB.iloc[0, 1] = "Y"
    aftab = pd.concat([aftab, aftabB])
    with pytest.raises(AssertionError) as e:
        HLAfreq.decrease_resolution(aftab, 2)
    assert "Multiple loci found for a single allele in a single population" in str(
        e.value
    )


def test_multi_sample_size_allele_decrease_resolution(aftab):
    aftabB = aftab.copy()
    aftabB.loc[0, "sample_size"] = 51
    aftab = pd.concat([aftab, aftabB])
    with pytest.raises(AssertionError) as e:
        HLAfreq.decrease_resolution(aftab, 2)
    assert (
        "Multiple sample_sizes found for a single allele in a single population"
        in str(e.value)
    )
