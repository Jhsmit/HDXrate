"""
This script is a python implementation of the Englander xls sheets for intrinsic rates
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

import numpy as np
from pathlib import Path

R = 1.987

# Activation energies (cal/mol)
E_act = {
    'acid': 14000.,
    'base': 17000.,
    'water': 19000.,
    'D': 1000.,
    'E': 1083.,
    'H': 7500.
}


def get_side_chain_dictionary(temperature, pH, k_reference):
    """
    Returns a dictionary with inductive effects of side chains on H/D exchange rates.

    Values are from [1]_, [2]_, [3]_, as described in [4]_

    Parameters
    ----------
    temperature: :obj:`float`
        Temperature in Kelvin.
    pH: :obj:`float`
        pH/pD value.
    k_reference: :obj:`dict`
        Dictionary of references values for amino acids D, E, H

    Returns
    -------

    constants: :obj:`dict`
        Dictionary of side chain modifiers. Values are: (acid_lambda, acid_rho, base_lambda, base_rho)

    References
    ----------

    .. [1] Bai, Y., Milne, J. S., Mayne, L. & Englander, S. W. Primary structure effects on peptide group hydrogen
       exchange. Proteins: Structure, Function, and Bioinformatics 17, 75–86 (1993).
    .. [2] Connelly, G. P., Bai, Y., Jeng, M.-F. & Englander, S. W. Isotope effects in peptide group hydrogen
       exchange. Proteins 17, 87–92 (1993).
    .. [3] Mori, S., Zijl, P. C. M. van & Shortle, D. Measurement of water–amide proton exchange rates in the
       denatured state of staphylococcal nuclease by a magnetization transfer technique. Proteins: Structure,
       Function, and Bioinformatics 28, 325–332 (1997).
    .. [4] Nguyen, D., Mayne, L., Phillips, M. C. & Walter Englander, S. Reference Parameters for Protein Hydrogen
       Exchange Rates. J. Am. Soc. Mass Spectrom. 29, 1936–1939 (2018).
    """

    root_dir = Path(__file__).parent
    names = ['name', 'short_name', 'acid_lambda', 'acid_rho', 'base_lambda', 'base_rho']
    side_chain_array = np.genfromtxt(root_dir / 'constants.txt', comments='#', skip_header=2, delimiter='\t', dtype=None,
                                     names=names, encoding=None, autostrip=True)

    side_chain_dict = {elem['short_name']: np.array(list(elem)[2:]) for elem in side_chain_array}
    for residue in ['D', 'E', 'H']: # residues D, E, H are calculated based on pH and pKa
        k_corrected = -np.log10(10**-k_reference[residue] * np.exp(-E_act[residue] * (1 / temperature - 1 / 278) / R)) # Check correct reference temperature

        deprotenated = side_chain_dict[residue + '0']
        protenated = side_chain_dict[residue + '+']

        values = np.log10(np.divide(10 ** (protenated - pH) + 10 ** (deprotenated - k_corrected),
                                    10 ** -k_corrected + 10 ** -pH))
        side_chain_dict[residue] = values
        if residue == 'E':
            side_chain_dict['CT'][0] = np.log10(np.divide(10 ** (0.05 - pH) + 10 ** (0.96 - k_corrected),
                                               10 ** -k_corrected + 10 ** -pH))

    return side_chain_dict


def correct_pH(pH_read, d_percentage=100.):
    """
    Correct for pH as described in Nguyen et al, 2018[1]_.
    This adds 0.4 to the pH value multiplied by the deuteration percentage.

    Note that there is no consensus on this correction factor. See also Rubinson, 2017[2]_

    Parameters
    ----------
    pH_read: :obj:`float`
        pH value of the solution as read by a standard glass electrode pH meter.
    d_percentage: :obj:`float`
        Percentage of deuterium in the solution.

    Returns
    -------
    pH_corrected : :obj:`float`
        Corrected pH value (pD)

    References
    ----------

   .. [1] Nguyen, D., Mayne, L., Phillips, M. C. & Walter Englander, S. Reference Parameters for Protein
      Hydrogen Exchange Rates. J. Am. Soc. Mass Spectrom. 29, 1936–1939 (2018).
   .. [2] Rubinson, K. A. Practical corrections for p(H,D) measurements in mixed H 2 O/D 2 O biological buffers.
      Anal. Methods 9, 2744–2750 (2017).
   """

    pH_corrected = pH_read + 0.4 * d_percentage / 100
    return pH_corrected


def k_int_from_sequence(sequence, temperature, pH_read, reference='poly', exchange_type='HD',
                        d_percentage=100., ph_correction=True, wildcard='X'):
    """
    Calculated intrisic rates of exchange for amide hydrogens in proteins.

    Calculations are based on  [1]_, [2]_, [3]_ and [4]_.


    Parameters
    ----------
    sequence: iterable object
        Input sequence in single-letter amino acid codes. Use 'Pc' for cis Proline, 'C2' for Cystine (disulfide)
    temperature: :obj:`float`
        Temperature of the exchange reaction in Kelvin.
    pH_read: :obj:`float`
        pH read by a standard glass electrode. If `ph_correction` is `True` this is corrected to pD in the case of `HD`
        exchange. pH changes due to temperature difference between measuring temperature and exchange temperature is
        buffer dependent and is not corrected for.
    reference: :obj:`str`
        Use `poly` or `oligo` reference data.
    exchange_type: :obj:`str`
        The type of exchange. Options are `HD`, `DH` or `HH`.
    d_percentage: :obj:`float`
        Percentage of Deuterium in the reaction solution. Used for pH/pD correction.
    ph_correction: :obj:`bool`
        Whether or not to correct the supplied `pH_read` value to pD. `DH` and `HH` exchange pH is not corrected.
    wildcard: :obj:`str`:
        Wildcard to use for unknown amino acids in the sequence. Amino acids with the wildcard and amino acids to the
        right (C-term) of wildcard residues return 0. as rate.

    Returns
    -------
    rates : :class:`~numpy.ndarray`
        Array with exchange rates in units of per second.

    References
    ----------

    .. [1] Bai, Y., Milne, J. S., Mayne, L. & Englander, S. W. Primary structure effects on peptide group hydrogen
       exchange. Proteins: Structure, Function, and Bioinformatics 17, 75–86 (1993).
    .. [2] Connelly, G. P., Bai, Y., Jeng, M.-F. & Englander, S. W. Isotope effects in peptide group hydrogen
       exchange. Proteins 17, 87–92 (1993).
    .. [3] Mori, S., Zijl, P. C. M. van & Shortle, D. Measurement of water–amide proton exchange rates in the
       denatured state of staphylococcal nuclease by a magnetization transfer technique. Proteins: Structure,
       Function, and Bioinformatics 28, 325–332 (1997).
    .. [4] Nguyen, D., Mayne, L., Phillips, M. C. & Walter Englander, S. Reference Parameters for Protein Hydrogen
       Exchange Rates. J. Am. Soc. Mass Spectrom. 29, 1936–1939 (2018).
    """

    if len(sequence) <3:
        raise ValueError('Sequence needs a minimum length of 3')
    if exchange_type not in ['HD', 'DH', 'HH']:
        raise ValueError(f"Unsupported exchange type '{exchange_type}'")

    if exchange_type == 'HD':
        exponents = np.array([1.62, 10.18, -1.5])
        pD = correct_pH(pH_read, d_percentage) if ph_correction else pH_read
        pKD = 15.05
        k_reference = {'D': 4.48, 'E': 4.93, 'H': 7.42}  # HD
    elif exchange_type == 'DH':
        exponents = np.array([1.4, 10., -1.6])
        pD = pH_read
        pKD = 14.17
        E_act['D'] -= 40
        k_reference = {'D': 3.87, 'E': 4.33, 'H': 7.0}  #DH
    elif exchange_type == 'HH':
        exponents = np.array([1.39, 10.08, -1.6])
        pD = pH_read
        pKD = 14.17
        E_act['D'] -= 40
        k_reference = {'D': 3.88, 'E': 4.35, 'H': 7.11}  #HH

    conc_D = 10. ** -pD
    conc_OD = 10. ** (pD - pKD)

    k_values = (10 ** exponents) / 60
    oligo_factors = [2.34, 1.35, 1.585]
    if reference == 'poly':
        k_acid_ref, k_base_ref, k_water_ref = k_values
    elif reference == 'oligo':
        k_acid_ref, k_base_ref, k_water_ref = k_values * oligo_factors
    else:
        raise ValueError("Value for 'reference' mush be either 'poly' or 'oligo'")

    sequence = list(sequence)
    sequence.insert(0, 'NT')
    sequence.append('CT')

    # Rates without inductive effects from neighbours, corrected for temperature
    k_acid = k_acid_ref * np.exp(-E_act['acid'] * (1 / temperature - 1 / 293) / R)
    k_base = k_base_ref * np.exp(-E_act['base'] * (1 / temperature - 1 / 293) / R)
    k_water = k_water_ref * np.exp(-E_act['water'] * (1 / temperature - 1 / 293) / R)

    side_chain_dict = get_side_chain_dictionary(temperature, pD, k_reference)

    k_int = []
    for i, residue in enumerate(sequence):
        if residue == 'NT':
            continue
        elif i == 1:  # First residue
            k_int.append(np.inf)
            continue

        next_residue = sequence[i + 1]
        prev_residue = sequence[i - 1]
        # Proline or unknown residues are set to zero rate
        if residue in ['P', 'Pc'] or wildcard in [prev_residue, residue]:
            k_int.append(0.)
            if next_residue == 'CT':
                break
            else:
                continue

        # Format is left_acid, right_acid, left_base, right_base
        _, prev_rho_acid, _, prev_rho_base = side_chain_dict[prev_residue]
        curr_lambda_acid, _, curr_lambda_base, _ = side_chain_dict[residue]

        if next_residue == 'CT':
            cterm_acid = side_chain_dict['CT'][0]
            cterm_base = side_chain_dict['CT'][2]

            Fa = 10 ** (prev_rho_acid + curr_lambda_acid + cterm_acid)
            Fb = 10 ** (curr_lambda_base + prev_rho_base + cterm_base)
        elif i == 2:  # Second residue in the chain (starts at 1)
            nterm_acid = side_chain_dict['NT'][1]
            nterm_base = side_chain_dict['NT'][3]

            Fa = 10 ** (prev_rho_acid + curr_lambda_acid + nterm_acid)
            Fb = 10 ** (curr_lambda_base + prev_rho_base + nterm_base)
        else:
            Fa = 10 ** (prev_rho_acid + curr_lambda_acid)
            Fb = 10 ** (curr_lambda_base + prev_rho_base)

        k_total_acid = Fa*k_acid*conc_D
        k_total_base = Fb*k_base*conc_OD
        k_total_water = Fb*k_water

        k_int.append(k_total_acid + k_total_base + k_total_water)

        if next_residue == 'CT':
            break

    return np.array(k_int)
