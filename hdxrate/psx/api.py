from hdxrate.psx.IntrinsicExchange import CalculateExchangeRatesForASingleChain
import numpy as np


def calc_k_int(sequence, temperature, pH, reference_data='poly'):
    """

    Parameters
    ----------
    sequence: iterable of strings
        Iterable of one-letter code amino acids. Unknown amino acids can be specified with 'X'
    temperature: :obj:`float`
        Temperature of the labelling buffer in Kelvin
    pH: :obj:`float`
        pH of labelling buffer. (pH read? pH corrected? pD?)
    reference_data: :obj:`str`
        'poly' or 'oligo'

    Returns
    -------

    k_int: ndarray, float
        Numpy array of intrinsic exchange rates (in units of per second)

    """
    sequence = list(sequence).copy()   # Function modifies list in-place
    k_int = CalculateExchangeRatesForASingleChain(sequence, temperature, pH, reference_data)


    return np.array(k_int)  # PSX returns intrinsic rates in per second


if __name__ == '__main__':
    # chain = list(
    #     'MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGTFPQLNLAPVNFDALFMNYLQQQAGEGTEEHQDA')

    chain = list('MSEQNNTEMTFQIQRIYTK')

    k_int = calc_k_int(chain, 300, 8., 'poly')
    np.array(k_int) * 60
    print(k_int)