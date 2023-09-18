#!/usr/bin/env python
import os
from rdkit import Chem
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

directory = "/Users/junfeng/datasets/CommercialLibraries/compoundlibrary_analysis"
sim_files = [
    "lifechemical_100k_simanalysis.sdf",
    "chemdiv_div_analysis.sdf",
    "wuxi_simanalysis.sdf"
    ]

colors = ["red","green","blue"]

tsne_files = [
    "lifechemical_skyhawk_tsne.sdf",
    "chemdiv_skyhawk_tsne.sdf",
    "wuxi_skyhawk_tsne.sdf"
  ]

plt.xlabel('Tanimoto')
plt.ylabel('Probability')
plt.title('Histogram of Similarity Distribution:')
for idx,f in enumerate(sim_files):
    tanimos = []
    #pic_name = os.path.splitext(f)[0]+".png"
    fname = os.path.join(directory,f)
    molfile = Chem.SDMolSupplier(fname)
    for mol in molfile:
        tanimoto = float(mol.GetProp("Tanimoto"))
        tanimos.append(tanimoto)
    sample = np.array(tanimos)
    stddev = np.std(sample)
    mean = np.mean(sample)
    weights = np.ones_like(sample) / (len(sample))
    # x = np.linspace(0, 1, 100)
    # plt.plot(x, norm.pdf(x, mean, stddev))
    n, bins, patches = plt.hist(sample, 20, weights=weights, facecolor=colors[idx], alpha=0.5)
    # y = norm.pdf(bins, mean, stddev)
    # plt.plot(bins, y, 'r--')
plt.savefig(os.path.join(directory, "test.png"),dpi=300)
