"""
Copyright 2019, University of Copenhagen

Martin Cramer Pedersen, mcpe@nbi.ku.dk

This file is part of PSX.

PSX is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License
as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

PSX is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with PSX. If not, see http://www.gnu.org/licenses/.
"""

###########
# Imports #
###########
import sys
import math


################################################################
# Function for determine exchange rates of an indivitual chain #
################################################################
def CalculateExchangeRatesForASingleChain(Chain, Temperature, pH, ReferenceData):
    # Get raference data
    if ReferenceData == "poly":
        ka = 10.0 ** 1.62
        kb = 10.0 ** 10.05
        kw = 10.0 ** -1.5

    elif ReferenceData == "oligo":
        ka = 10.0 ** 1.62 * 2.34
        kb = 10.0 ** 10.05 * 1.35
        kw = 10.0 ** -1.5 * 1.585

    else:
        sys.stderr.write("Either poly or oligo must be selected as reference data. Exiting.\n")
        sys.exit()

    # Temperature correction
    R = 1.987
    TemperatureCorrection = (1.0 / Temperature - 1.0 / 293.0) / R

    # Activation energies (in cal/mol)
    AcidActivationEnergy = 14000.0
    BaseActivationEnergy = 17000.0
    SolventActivationEnergy = 19000.0

    AspActivationEnergy = 1000.0
    GluActivationEnergy = 1083.0
    HisActivationEnergy = 7500.0

    # Corrections based on activation energies
    AcidTemperatureCorrection = math.exp(- TemperatureCorrection * AcidActivationEnergy)
    BaseTemperatureCorrection = math.exp(- TemperatureCorrection * BaseActivationEnergy)
    SolventTemperatureCorrection = math.exp(- TemperatureCorrection * SolventActivationEnergy)

    AspTemperatureCorrection = math.exp(- TemperatureCorrection * AspActivationEnergy)
    GluTemperatureCorrection = math.exp(- TemperatureCorrection * GluActivationEnergy)
    HisTemperatureCorrection = math.exp(- TemperatureCorrection * HisActivationEnergy)

    # Corrected pH in D2O
    pH += 0.4

    # pK-values
    pKD = 15.05
    pKAsp = 4.48 * AspTemperatureCorrection
    pKGlu = 4.93 * GluTemperatureCorrection
    pKHis = 7.40 * HisTemperatureCorrection

    LambdaProtonatedAcidAsp = math.log(
        10.0 ** (-0.90 - pH) / (10.0 ** -pKAsp + 10.0 ** -pH) + 10.0 ** (0.90 - pKAsp) / (10.0 ** -pKAsp + 10.0 ** -pH),
        10.0)
    LambdaProtonatedAcidGlu = math.log(
        10.0 ** (-0.60 - pH) / (10.0 ** -pKGlu + 10.0 ** -pH) + 10.0 ** (-0.90 - pKGlu) / (
                10.0 ** -pKGlu + 10.0 ** -pH), 10.0)
    LambdaProtonatedAcidHis = math.log(
        10.0 ** (-0.80 - pH) / (10.0 ** -pKHis + 10.0 ** -pH) + 10.0 ** (0.00 - pKHis) / (10.0 ** -pKHis + 10.0 ** -pH),
        10.0)

    RhoProtonatedAcidAsp = math.log(
        10.0 ** (-0.12 - pH) / (10.0 ** -pKAsp + 10.0 ** -pH) + 10.0 ** (0.58 - pKAsp) / (10.0 ** -pKAsp + 10.0 ** -pH),
        10.0)
    RhoProtonatedAcidGlu = math.log(
        10.0 ** (-0.27 - pH) / (10.0 ** -pKGlu + 10.0 ** -pH) + 10.0 ** (0.31 - pKGlu) / (10.0 ** -pKGlu + 10.0 ** -pH),
        10.0)
    RhoProtonatedAcidHis = math.log(
        10.0 ** (-0.51 - pH) / (10.0 ** -pKHis + 10.0 ** -pH) + 10.0 ** (0.00 - pKHis) / (10.0 ** -pKHis + 10.0 ** -pH),
        10.0)

    LambdaProtonatedBaseAsp = math.log(
        10.0 ** (0.69 - pH) / (10.0 ** -pKAsp + 10.0 ** -pH) + 10.0 ** (0.10 - pKAsp) / (10.0 ** -pKAsp + 10.0 ** -pH),
        10.0)
    LambdaProtonatedBaseGlu = math.log(
        10.0 ** (0.24 - pH) / (10.0 ** -pKGlu + 10.0 ** -pH) + 10.0 ** (-0.11 - pKGlu) / (10.0 ** -pKGlu + 10.0 ** -pH),
        10.0)
    LambdaProtonatedBaseHis = math.log(
        10.0 ** (0.80 - pH) / (10.0 ** -pKHis + 10.0 ** -pH) + 10.0 ** (-0.10 - pKHis) / (10.0 ** -pKHis + 10.0 ** -pH),
        10.0)

    RhoProtonatedBaseAsp = math.log(
        10.0 ** (0.60 - pH) / (10.0 ** -pKAsp + 10.0 ** -pH) + 10.0 ** (-0.18 - pKAsp) / (10.0 ** -pKAsp + 10.0 ** -pH),
        10.0)
    RhoProtonatedBaseGlu = math.log(
        10.0 ** (0.39 - pH) / (10.0 ** -pKGlu + 10.0 ** -pH) + 10.0 ** (-0.15 - pKGlu) / (10.0 ** -pKGlu + 10.0 ** -pH),
        10.0)
    RhoProtonatedBaseHis = math.log(
        10.0 ** (0.83 - pH) / (10.0 ** -pKHis + 10.0 ** -pH) + 10.0 ** (0.14 - pKHis) / (10.0 ** -pKHis + 10.0 ** -pH),
        10.0)

    # Termini
    RhoAcidNTerm = -1.32
    LambdaAcidCTerm = math.log(
        10.0 ** (0.05 - pH) / (10.0 ** -pKGlu + 10.0 ** -pH) + 10.0 ** (0.96 - pKGlu) / (10.0 ** -pKGlu + 10.0 ** -pH),
        10.0)

    RhoBaseNTerm = 1.62
    LambdaBaseCTerm = -1.80

    # Ion concentrations
    DIonConc = 10.0 ** -pH
    ODIonConc = 10.0 ** (pH - pKD)

    # Dictionary for acid values (format is (lambda, rho))
    MilneAcid = {}

    MilneAcid["NTerminal"] = (None, RhoAcidNTerm)
    MilneAcid["CTerminal"] = (LambdaAcidCTerm, None)

    MilneAcid["A"] = (0.00, 0.00)
    MilneAcid["C"] = (-0.54, -0.46)
    MilneAcid["C2"] = (-0.74, -0.58)
    MilneAcid["D"] = (LambdaProtonatedAcidAsp, RhoProtonatedAcidAsp)
    MilneAcid["D+"] = (-0.90, -0.12)
    MilneAcid["E"] = (LambdaProtonatedAcidGlu, RhoProtonatedAcidGlu)
    MilneAcid["E+"] = (-0.60, -0.27)
    MilneAcid["F"] = (-0.52, -0.43)
    MilneAcid["G"] = (-0.22, 0.22)
    MilneAcid["H"] = (LambdaProtonatedAcidHis, RhoProtonatedAcidHis)
    MilneAcid["I"] = (-0.91, -0.59)
    MilneAcid["K"] = (-0.56, -0.29)
    MilneAcid["L"] = (-0.57, -0.13)
    MilneAcid["M"] = (-0.64, -0.28)
    MilneAcid["N"] = (-0.58, -0.13)
    MilneAcid["P"] = (0.00, -0.19)
    MilneAcid["Pc"] = (0.00, -0.85)
    MilneAcid["Q"] = (-0.47, -0.27)
    MilneAcid["R"] = (-0.59, -0.32)
    MilneAcid["S"] = (-0.44, -0.39)
    MilneAcid["T"] = (-0.79, -0.47)
    MilneAcid["V"] = (-0.74, -0.30)
    MilneAcid["W"] = (-0.40, -0.44)
    MilneAcid["Y"] = (-0.41, -0.37)

    # Dictionary for base values (format is (lambda, rho))
    MilneBase = {}

    MilneBase["NTerminal"] = (None, RhoBaseNTerm)
    MilneBase["CTerminal"] = (LambdaBaseCTerm, None)

    MilneBase["A"] = (0.00, 0.00)
    MilneBase["C"] = (0.62, 0.55)
    MilneBase["C2"] = (0.55, 0.46)
    MilneBase["D"] = (LambdaProtonatedBaseAsp, RhoProtonatedBaseAsp)
    MilneBase["D+"] = (0.69, 0.60)
    MilneBase["E"] = (LambdaProtonatedBaseGlu, RhoProtonatedBaseGlu)
    MilneBase["E+"] = (0.24, 0.39)
    MilneBase["F"] = (-0.24, 0.06)
    MilneBase["G"] = (0.27, 0.17)
    MilneBase["H"] = (LambdaProtonatedBaseHis, RhoProtonatedBaseHis)
    MilneBase["I"] = (-0.73, -0.23)
    MilneBase["K"] = (-0.04, 0.12)
    MilneBase["L"] = (-0.58, -0.21)
    MilneBase["M"] = (-0.01, 0.11)
    MilneBase["N"] = (0.49, 0.32)
    MilneBase["P"] = (0.00, -0.24)
    MilneBase["Pc"] = (0.00, 0.60)
    MilneBase["Q"] = (0.06, 0.20)
    MilneBase["R"] = (0.08, 0.22)
    MilneBase["S"] = (0.37, 0.30)
    MilneBase["T"] = (-0.07, 0.20)
    MilneBase["V"] = (-0.70, -0.14)
    MilneBase["W"] = (-0.41, -0.11)
    MilneBase["Y"] = (-0.27, 0.05)

    # Default values
    MilneAcid["?"] = (0.00, 0.00)
    MilneBase["?"] = (0.00, 0.00)

    # Loop over the chain
    IntrinsicEnchangeRates = [0.0]
    Chain.insert(0, "NTerminal")
    Chain.append("CTerminal")

    # Account for middle residues
    # Iteration starts at two, +1 for appending of "NTerminal", +1 for skipping the N terminal residue
    for i in range(2, len(Chain) - 1):
        Residue = Chain[i]

        LeftResidue = Chain[i - 1]
        RightResidue = Chain[i + 1]

        if Residue in ("P", "Pc"):
            IntrinsicEnchangeRates.append(0.0)
        # If the residue or the residue on the left is unknown, skip it and add 0.0
        elif 'X' in [Residue, LeftResidue]:  # 'X' added by JHS
            IntrinsicEnchangeRates.append(0.0)
        else:
            # Identify neighbors

            if RightResidue == "CTerminal":
                Fa = 10.0 ** (MilneAcid[LeftResidue][1] + MilneAcid[Residue][0] + MilneAcid["CTerminal"][0])
                Fb = 10.0 ** (MilneBase[LeftResidue][1] + MilneBase[Residue][0] + MilneBase["CTerminal"][0])

            elif i == 2:  # Second residue in the chain
                Fa = 10.0 ** (MilneAcid["NTerminal"][1] + MilneAcid[LeftResidue][1] + MilneAcid[Residue][0])
                Fb = 10.0 ** (MilneBase["NTerminal"][1] + MilneBase[LeftResidue][1] + MilneBase[Residue][0])

            else:
                Fa = 10.0 ** (MilneAcid[LeftResidue][1] + MilneAcid[Residue][0])
                Fb = 10.0 ** (MilneBase[LeftResidue][1] + MilneBase[Residue][0])

            # Contributions from acid, base, and water
            kaT = Fa * ka * AcidTemperatureCorrection * DIonConc
            kbT = Fb * kb * BaseTemperatureCorrection * ODIonConc
            kwT = Fb * kw * SolventTemperatureCorrection

            # Collect exchange rates
            IntrinsicExchangeRate = kaT + kbT + kwT

            # Construct list
            IntrinsicEnchangeRates.append(IntrinsicExchangeRate)

    # Rescale from 1/min to 1/s
    IntrinsicEnchangeRates = [IntrinsicEnchangeRates[i] / 60.0 for i in range(len(IntrinsicEnchangeRates))]

    return IntrinsicEnchangeRates

