from hdxrate.expfact.kint import calculate_kint_per_residue
import numpy as np


def calc_k_int(sequence, temperature, pH):
    """

    Parameters
    ----------
    sequence: iterable of strings
        Iterable of one-letter code amino acids. Unknown amino acids can be specified with 'X'
    temperature: :obj:`float`
        Temperature of the labelling buffer in Kelvin
    pH: :obj:`float`
        pH of labelling buffer. (pH read? pH corrected? pD?)

    Returns
    -------

    k_int: ndarray, float
        Numpy array of intrinsic exchange rates (in units of per second)

    """

    list(sequence).copy()

    k_int_list = [-1.]  # first residue
    for i, (previous, current) in enumerate(zip(sequence[:-1], sequence[1:])):
        if previous == 'X' or current == 'X':
            k_int_list.append(0.)
        elif current == 'P':
            k_int_list.append(0.)
        else:
            k_int = calculate_kint_per_residue(previous, current, i + 2, len(sequence), temperature, pH)
            k_int_list.append(k_int)

    return np.array(k_int_list) / 60


if __name__ == '__main__':
    chain = list('MSEQNNTEMTFQIQRIYTK')
    k_int = calc_k_int(chain, 300, 8.)

    print(k_int)
