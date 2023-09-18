#!/usr/bin/env python
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import sys
#from tsnecuda import TSNE
import umap,hdbscan

def tanimoto_dist(a,b):
    dotprod = np.dot(a,b)
    tc = dotprod / (np.sum(a) + np.sum(b) - dotprod)
    return 1.0-tc

def calc_fp_arr( mols ):
    fplist = []
    for mol in mols:
        arr = np.zeros((1,))
        fp = AllChem.GetMorganFingerprint(mol,3)
        #fp = AllChem.GetAtomPairFingerprint(mol)
        DataStructs.ConvertToNumpyArray(fp, arr)
        fplist.append(arr)
    return np.asarray(fplist)


sdf = sys.argv[1]
drugs = [ mol for mol in Chem.SDMolSupplier(sdf) if mol != None ]

res =calc_fp_arr(drugs)

# umap.UMAP(
#     n_neighbors=30,
#     min_dist=0.0,
#     n_components=2,
#     random_state=42,
# ).fit_transform(mnist.data)

clustering = umap.UMAP(metric="jaccard",n_neighbors=7,n_components=3).fit_transform(res)
print(clustering.shape)
X = clustering[:, 0]
Y = clustering[:, 1]
#Z = clustering[:, 2]
#labels = hdbscan.HDBSCAN().fit_predict(clustering)
labels = DBSCAN().fit_predict(clustering)

sd_writer = Chem.SDWriter(sys.argv[2])
n_unclustered = 0
for id,mol in enumerate(drugs):
    if labels[id] == -1:
        n_unclustered += 1
    mol.SetProp("Cluster",str(labels[id]))
    mol.SetProp("X_umap",str(X[id]))
    mol.SetProp("Y_umap",str(Y[id]))
#   mol.SetProp("Z_umap",str(Z[id]))
    sd_writer.write(mol)
sd_writer.close()
print("UnClustered compounds:%d"%n_unclustered)
print(len(np.unique(labels)))

