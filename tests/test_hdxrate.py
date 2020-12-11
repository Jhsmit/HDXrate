"""Tests for `hdxrate` package."""

import numpy as np
from hdxrate import k_int_from_sequence


class TestHDXrate(object):
    """Tests for `hdxrate` package."""

    @classmethod
    def setup_class(cls):
        cls.sequence = list('AAAWADEAA')

    def test_k_int_HD(self):
        pH_read = 6.6
        temp = 279

        rates = k_int_from_sequence(self.sequence, temp, pH_read, reference='poly',
                                    ph_correction='englander', wildcard='X') * 60
        # Reference rates obtained from Englander group xls sheet
        ref_rates = np.array([np.inf, 1.29939811E+03, 3.11703908E+01, 1.21266892E+01, 2.41959255E+01, 3.95805093E+01,
                              1.63948783E+01, 2.25232682E+01, 4.94028674E-01])
        assert np.allclose(rates, ref_rates)

    def test_k_int_DH(self):
        pH_read = 7.0
        temp = 278

        rates = k_int_from_sequence(self.sequence, temp, pH_read, exchange_type='DH', reference='poly',
                                    ph_correction='englander', wildcard='X') * 60
        # Reference rates obtained from Englander group xls sheet
        ref_rates = np.array([np.inf, 5.831286065E+03, 1.398828104E+02, 5.442072826E+01, 1.085836280E+02,
                              1.764790873E+02, 7.219785837E+01, 9.955072660E+01, 2.216999376E+00])
        assert np.allclose(rates, ref_rates)

    def test_k_int_HH(self):
        pH_read = 7.0
        temp = 278

        rates = k_int_from_sequence(self.sequence, temp, pH_read, exchange_type='HH', reference='poly',
                                    ph_correction='englander', wildcard='X') * 60
        # Reference rates obtained from Englander group xls sheet
        ref_rates = np.array([np.inf, 7.01071144E+03, 1.68175255E+02, 6.54277663E+01, 1.30545556E+02, 2.12183978E+02,
                              8.68187188E+01, 1.19715146E+02, 2.66540425E+00])
        assert np.allclose(rates, ref_rates)
