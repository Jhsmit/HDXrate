"""
This files calculates and plots intrinsic rates of randomly generated protein sequences and looks at differences between
methods.

"""

from hdxrate.expfact.api import calc_k_int as expfact_k_int
from hdxrate.psx.api import calc_k_int as psx_k_int
import numpy as np
import matplotlib.pyplot as plt

amino_acids = list('GALMFWKQESPVICYHRNDT')
repeats = 100
temperature = 300
pHs = [6, 7, 8, 9]
fig, axes = plt.subplots(2, 2)

for ax, pH in zip(axes.flatten(), pHs):
    diffs = []
    for i in range(repeats):
        sequence = list(np.random.choice(np.array(amino_acids), 800))
        r_number = np.arange(len(sequence)) + 1

        psx_poly = psx_k_int(sequence, temperature, pH, 'poly')
        expfact = expfact_k_int(sequence, temperature, pH)

        d = psx_poly / expfact
        diffs += list(psx_poly / expfact)

    diffs = np.array(diffs)
    diffs = diffs[~np.isnan(diffs)]

    ax.set_yscale('log')
    ax.hist(diffs, bins=75)
    ax.set_title(f'pH {pH}')
    ax.set_ylabel('Count')
    ax.set_xlabel('Fold difference')
#fig.suptitle('Intrisic rate differences expfact / psx')
plt.tight_layout()
plt.savefig('Rate differences histograms.png')


# #
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 5), sharex=True)
# ax1.set_yscale('log')
# ax1.scatter(r_number, psx_poly, label='psx_poly')
# ax1.scatter(r_number, psx_oligo, label='psx_oligo')
# ax1.scatter(r_number, expfact, label='expfact')
# ax1.set_ylim(0.1, 1e5)
# ax1.set_ylabel('Intrinsic Exchange \n rate (1/s)')
# ax1.legend()
#
# ax2.plot(r_number, psx_poly / expfact, color='k', label='psx_poly / expfact')
# ax2.plot(r_number, psx_oligo / expfact, color='k', linestyle='--', label='psx_oligo / expfact')
# ax2.set_ylabel('Fold difference')
# ax2.set_xlabel('Residue number')
# ax2.legend()
# ax2.plot()
# plt.show()
#
# #
# # fig, ax = plt.subplots()
# # #ax.set_yscale('log')
# # #ax.plot(psx, label='psx')
# # ax.plot(expfact / psx, label='expfact')
# # ax.legend()
# # plt.show()