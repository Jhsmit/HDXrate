"""Tests for `hdxrate` package."""

import numpy as np
from hdxrate import k_int_from_sequence


class TestHDXrate(object):
    """Tests for `hdxrate` package."""

    @classmethod
    def setup_class(cls):
        cls.sequence = list('AAAWADEAA')

    def test_k_int_calculation(self):
        seq = 'AAAWADEAA'

        pH_read = 6.6
        temp = 279

        rates = k_int_from_sequence(seq, temp, pH_read, reference='poly', ph_correction='englander', wildcard='X') * 60
        ref_rates = np.array([np.inf, 1.29939811E+03, 3.11703908E+01, 1.21266892E+01, 2.41959255E+01, 3.95805093E+01,
                              1.63948783E+01, 2.25232682E+01, 4.94028674E-01])

        assert np.allclose(rates[1:], ref_rates[1:])

