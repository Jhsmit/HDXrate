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

import numpy as np
from pathlib import Path

R = 1.987
pKD = 15.05

# Activation energies (cal/mol)
E_act = {
    'acid': 14000.,
    'base': 17000.,
    'water': 19000.,
    'D': 1000.,
    'E': 1083.,
    'H': 7500.
}

root_dir = Path(__file__).parent

def get_side_chain_dictionary(temperature, pH):
    """

    Parameters
    ----------
    temperature
    pH

    Returns
    -------

    constants: :obj:`dict`
        List of side chain modifiers (acid_lambda, acid_rho, base_lambda, base_rho)

    """

    names = ['name', 'short_name', 'acid_lambda', 'acid_rho', 'base_lambda', 'base_rho']
    side_chain_array = np.genfromtxt(root_dir / 'constants.txt', comments='#', skip_header=2, delimiter='\t', dtype=None,
                                     names=names, encoding=None, autostrip=True)

    side_chain_dict = {elem['short_name']: np.array(list(elem)[2:]) for elem in side_chain_array}
    for residue in ['D', 'E', 'H']: # residues D, E, H are calculated based on pH and pKa
        k_reference = {'D': 4.48, 'E': 4.93, 'H': 7.42}
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


def correct_pH(pH_read, method='englander'):
    return pH_read - 0.4


def k_int_from_sequence(sequence, temperature, pH_read, reference='poly', ph_correction='englander', wildcard='X'):
    # Reference rates Ala-Ala corrected to 20C in units of per second.
    # Nguyen et al, 2018, Table 1. / Englander xls sheets
    #todo update for DH exchange

    if len(sequence) <3:
        raise ValueError('Sequence needs a minimum length of 3')

    pD = correct_pH(pH_read, method=ph_correction)
    conc_D = 10.**-pD
    conc_OD = 10.**(pD - pKD)

    exponents = np.array([1.62, 10.18, -1.5])
    k_values = (10 ** exponents) / 60

    if reference == 'poly':
        k_acid_ref, k_base_ref, k_water_ref = k_values
    elif reference == 'oligo':
        oligo_factors = [2.34, 1.35, 1.585]
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

    side_chain_dict = get_side_chain_dictionary(temperature, pD)

    k_int = []
    for i, residue in enumerate(sequence):
        if residue == 'NT':
            continue
        elif i == 1:
            k_int.append(np.inf)
            continue

        next_residue = sequence[i + 1]
        prev_residue = sequence[i - 1]
        if residue in ['P', 'Pc']:  # Proline residues do not exchange
            k_int.append(0.)
            continue
        if wildcard in [prev_residue, residue]:  # Unknown residues get assigned rate 0
            k_int.append(0.)
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
