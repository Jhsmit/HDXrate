"""Main module."""

from hdxrate.expfact.api import calc_k_int as expfact_k_int
from hdxrate.psx.api import calc_k_int as psx_k_int


def k_int_from_sequence(sequence, temperature, pH, module='expfact', **kwargs):
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
        Which module to use for calculating intrinsic rates. Default is 'expfact', options are 'expfact', 'psx'
    **kwargs:
        Additional module-specific kwargs to pass on intrinsic rate calculation. See module .api for details.

    Returns
    -------

    k_int: ndarray, float
        Numpy array of intrinsic exchange rates (in units of per second)

    """

    sequence = list(sequence)
    if module == 'expfact':
        func = expfact_k_int
    elif module == 'psx':
        func = psx_k_int
    else:
        raise ValueError(f"Invalid value '{module}' specified for module. Options are 'expfact' or 'psx')")

    k_int = func(sequence, temperature, pH, **kwargs)

    return k_int