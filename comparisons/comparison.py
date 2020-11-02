"""
This files calculates and plots intrinsic rates of a sample protein (ecSecB) and looks at differences between
methods.
"""

from hdxrate.expfact.api import calc_k_int as expfact_k_int
from hdxrate.psx.api import calc_k_int as psx_k_int
import numpy as np
import matplotlib.pyplot as plt


sequence = list('MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLRVTVTASLG')
r_number = np.arange(len(sequence)) + 1

temperature = 300
pH = 8.

psx_poly = psx_k_int(sequence, temperature, pH, 'poly')
psx_oligo = psx_k_int(sequence, temperature, pH, 'oligo')
expfact = expfact_k_int(sequence, temperature, pH)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 5), sharex=True)
ax1.set_yscale('log')
ax1.plot(r_number, psx_poly, label='psx_poly', marker='o', linestyle='None')
ax1.plot(r_number, psx_oligo, label='psx_oligo', marker='o', linestyle='None')
ax1.plot(r_number, expfact, label='expfact', marker='o', linestyle='None')
#ax1.set_ylim(0.1, 1e5)
ax1.set_xlim(0, 80)
ax1.set_ylabel('Intrinsic Exchange \n rate (1/s)')
ax1.legend()

ax2.plot(r_number, psx_poly / expfact, color='k', label='psx_poly / expfact')
ax2.plot(r_number, psx_oligo / expfact, color='k', linestyle='--', label='psx_oligo / expfact')
ax2.set_ylabel('Fold difference')
ax2.set_xlabel('Residue number')
ax2.legend()
plt.savefig('Rate differences.png')