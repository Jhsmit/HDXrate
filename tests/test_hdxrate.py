"""Tests for `hdxrate` package."""

import numpy as np
from hdxrate import k_int_from_sequence
from hdxrate.hdxrate import get_side_chain_dictionary
from pathlib import Path
from functools import reduce
from itertools import combinations
from operator import add

pth = Path(__file__).parent

class TestHDXrate(object):
    """Tests for `hdxrate` package."""

    @classmethod
    def setup_class(cls):
        cls.seq1 = list('AAAWADEAA')

        k_reference = {'D': 3.87, 'E': 4.33, 'H': 7.0}  # DH
        dict = get_side_chain_dictionary(278, 8, k_reference)
        one_letter = [k for k in dict.keys() if len(k) == 1]
        cls.seq2 = reduce(add, [a + b for a, b in combinations(one_letter, 2)])

        cls.seq2 = np.genfromtxt(pth / 'sequence.txt', dtype='U')

    def test_seq1(self):
        # HD rates
        pH_read = 6.6
        temp = 279
        rates = k_int_from_sequence(self.seq1, temp, pH_read, reference='poly', wildcard='X') * 60
        # Reference rates obtained from Englander group xls sheet
        ref_rates = np.array([np.inf, 1.29939811E+03, 3.11703908E+01, 1.21266892E+01, 2.41959255E+01, 3.95805093E+01,
                              1.63948783E+01, 2.25232682E+01, 4.94028674E-01])
        assert np.allclose(rates, ref_rates)

        # DH rates
        pH_read = 7.0
        temp = 278
        rates = k_int_from_sequence(self.seq1, temp, pH_read, exchange_type='DH', reference='poly') * 60
        # Reference rates obtained from Englander group xls sheet
        ref_rates = np.array([np.inf, 5.831286065E+03, 1.398828104E+02, 5.442072826E+01, 1.085836280E+02,
                              1.764790873E+02, 7.219785837E+01, 9.955072660E+01, 2.216999376E+00])
        assert np.allclose(rates, ref_rates)

        #HH rates
        pH_read = 7.0
        temp = 278
        rates = k_int_from_sequence(self.seq1, temp, pH_read, exchange_type='HH', reference='poly', wildcard='X') * 60
        # Reference rates obtained from Englander group xls sheet
        ref_rates = np.array([np.inf, 7.01071144E+03, 1.68175255E+02, 6.54277663E+01, 1.30545556E+02, 2.12183978E+02,
                              8.68187188E+01, 1.19715146E+02, 2.66540425E+00])
        assert np.allclose(rates, ref_rates)

    def test_seq2(self):
        reference_rates = np.genfromtxt(pth / 'exchange_rates_xls.txt', skip_header=2,
                                        delimiter='\t', filling_values=0.)

        reference_rates[0][:] = np.inf

        #HD exchange
        rates = k_int_from_sequence(self.seq2, 279, 7 - 0.4, exchange_type='HD') * 60
        assert np.allclose(rates, reference_rates[:, 0], rtol=0.1, equal_nan=True)

        #
        rates = k_int_from_sequence(self.seq2, 279, 7, exchange_type='DH') * 60
        assert np.allclose(rates, reference_rates[:, 1], rtol=0.1, equal_nan=True)

        rates = k_int_from_sequence(self.seq2, 279, 7, exchange_type='HH') * 60
        assert np.allclose(rates, reference_rates[:, 2], rtol=0.1, equal_nan=True)