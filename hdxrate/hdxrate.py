"""
This script provides and API to psx intrinsic exchange rate calculation
Copyright (C) 2020 Jochem Smit

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""


from hdxrate.psx.api import calc_k_int as psx_k_int


def k_int_from_sequence(sequence, temperature, pH, module='psx', **kwargs):
    """
    Calculate intrinsic rate of amide hydrogen exchange.

    The first and last entries in the 'sequence' iterable are assumed to be the N- and C-terminal residues, respectively,
    which affects their rate of exchange. Sequences can be padded with 'X' entries to correct this.


    Parameters
    ----------
    sequence: iterable of strings
        Iterable of one-letter code amino acids. Unknown amino acids can be specified with 'X'
    temperature: :obj:`float`
        Temperature of the labelling buffer in Kelvin
    pH: :obj:`float`
        pH of labelling buffer. (pH read? pH corrected? pD?)
    module: :obj:`str`
        Which module to use for calculating intrinsic rates. Default is 'psx'.
    **kwargs:
        Additional module-specific kwargs to pass on intrinsic rate calculation. See module .api for details.

    Returns
    -------

    k_int: ndarray, float
        Numpy array of intrinsic exchange rates (in units of per second)

    """

    sequence = list(sequence)
    if module == 'psx':
        func = psx_k_int
    else:
        raise ValueError(f"Invalid value '{module}' specified for module. Options are 'psx'")

    k_int = func(sequence, temperature, pH, **kwargs)

    return k_int