"""Tests for `hdxrate` package."""

import numpy as np
from hdxrate.expfact.api import calc_k_int as expfact_k_int
from hdxrate.psx.api import calc_k_int as psx_k_int
from hdxrate import k_int_from_sequence


class TestHdxrate(object):
    """Tests for `hdxrate` package."""

    @classmethod
    def setup_class(cls):
        cls.sequence = list('MSEQNNTEMTFQIQRIYTK')

    def test_k_int_calculation(self):

        expfact = expfact_k_int(self.sequence, 300, 8.)
        k_int = k_int_from_sequence(self.sequence, 300, 8., 'expfact')
        assert np.allclose(k_int, expfact)

        psx = psx_k_int(self.sequence, 300, 8., 'poly')
        k_int = k_int_from_sequence(self.sequence, 300, 8., 'psx')

        assert np.allclose(k_int, psx)

