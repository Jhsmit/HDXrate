"""Tests for `hdxrate` package."""

import numpy as np
from hdxrate import k_int_from_sequence
from hdxrate.hdxrate import get_side_chain_dictionary
from pathlib import Path
from functools import reduce
from itertools import combinations
from operator import add

import pytest

pth = Path(__file__).parent


@pytest.fixture()
def seq1():
    return list("AAAWADEAA")


@pytest.fixture()
def seq2():
    """sequence two a sequence of the pairwise combination of all side chains"""
    k_reference = {"D": 3.87, "E": 4.33, "H": 7.0}  # DH
    chains_dict = get_side_chain_dictionary(278, 8, k_reference)
    one_letter = [k for k in chains_dict.keys() if len(k) == 1]
    seq2 = reduce(add, [a + b for a, b in combinations(one_letter, 2)])

    return list(seq2)


def test_seq1_HD(seq1):
    # HD rates
    pH_read = 6.6
    temp = 279
    rates = (
        k_int_from_sequence(seq1, temp, pH_read, reference="poly", wildcard="X") * 60
    )
    # Reference rates obtained from Englander group xls sheet
    ref_rates = np.array(
        [
            np.inf,
            1.29939811e03,
            3.11703908e01,
            1.21266892e01,
            2.41959255e01,
            3.95805093e01,
            1.63948783e01,
            2.25232682e01,
            4.94028674e-01,
        ]
    )
    assert np.allclose(rates, ref_rates)


def test_seq1_DH(seq1):
    # DH rates
    pH_read = 7.0
    temp = 278
    rates = (
        k_int_from_sequence(seq1, temp, pH_read, exchange_type="DH", reference="poly")
        * 60
    )
    # check repeated calls give the same result
    rep_rates = (
        k_int_from_sequence(seq1, temp, pH_read, exchange_type="DH", reference="poly")
        * 60
    )
    # Reference rates obtained from Englander group xls sheet
    ref_rates = np.array(
        [
            np.inf,
            5.831286065e03,
            1.398828104e02,
            5.442072826e01,
            1.085836280e02,
            1.764790873e02,
            7.219785837e01,
            9.955072660e01,
            2.216999376e00,
        ]
    )
    assert np.allclose(rates, ref_rates)
    assert np.allclose(rates, rep_rates)


def test_seq1_HH(seq1):
    # HH rates
    pH_read = 7.0
    temp = 278
    rates = (
        k_int_from_sequence(
            seq1,
            temp,
            pH_read,
            exchange_type="HH",
            reference="poly",
            wildcard="X",
        )
        * 60
    )

    rep_rates = (
        k_int_from_sequence(
            seq1,
            temp,
            pH_read,
            exchange_type="HH",
            reference="poly",
            wildcard="X",
        )
        * 60
    )

    # Reference rates obtained from Englander group xls sheet
    ref_rates = np.array(
        [
            np.inf,
            7.01071144e03,
            1.68175255e02,
            6.54277663e01,
            1.30545556e02,
            2.12183978e02,
            8.68187188e01,
            1.19715146e02,
            2.66540425e00,
        ]
    )
    assert np.allclose(rates, ref_rates)
    assert np.allclose(rates, rep_rates)


def test_seq2(seq2):
    reference_rates = np.genfromtxt(
        pth / "exchange_rates_xls.txt",
        skip_header=2,
        delimiter="\t",
        filling_values=0.0,
    )

    reference_rates[0][:] = np.inf

    # HD exchange
    rates = k_int_from_sequence(seq2, 279, 7 - 0.4, exchange_type="HD") * 60
    assert np.allclose(rates, reference_rates[:, 0], rtol=0.1, equal_nan=True)

    # DH exchange
    rates = k_int_from_sequence(seq2, 279, 7, exchange_type="DH") * 60
    assert np.allclose(rates, reference_rates[:, 1], rtol=0.1, equal_nan=True)

    # HH exchange
    rates = k_int_from_sequence(seq2, 279, 7, exchange_type="HH") * 60
    assert np.allclose(rates, reference_rates[:, 2], rtol=0.1, equal_nan=True)
