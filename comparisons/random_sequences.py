"""
This files calculates and plots intrinsic rates of randomly generated protein sequences and looks at differences between
methods.

"""

from hdxrate.expfact.api import calc_k_int as expfact_k_int
from hdxrate.psx.api import calc_k_int as psx_k_int
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(43)
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
    diffs = diffs[~np.isnan(diffs)]  # remove NaN

    ax.set_yscale('log')
    ax.hist(diffs, bins=75)
    ax.set_title(f'pH {pH}')
    ax.set_ylabel('Count')
    ax.set_xlabel('Fold difference')
    ax.set_xlim(0, None)
    #print(pH, np.min(diffs), np.max(diffs)) :
    # 6 1.586076428945057 4.323111639721702
    # 7 1.6062294252044966 6.2163497448764105
    # 8 1.606402448281534 31.978492966181015
    # 9 1.6064194484834653 45.168683743072926

plt.tight_layout()
#plt.show()
plt.savefig('Rate differences histograms.png')
