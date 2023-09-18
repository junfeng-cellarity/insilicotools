from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs

import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import sys

if len(sys.argv) != 4:
    print("Usage:%s input.sdf output.sdf numClusters"%sys.argv[0])
    sys.exit(1)

drugs = [ mol for mol in Chem.SDMolSupplier( sys.argv[1] ) if mol != None ]
nClusters = int(sys.argv[3])

def calc_fp_arr( mols ):
    fplist = []
    for mol in mols:
        arr = np.zeros( (1,) )
        fp = AllChem.GetMorganFingerprintAsBitVect( mol, 3 )
        DataStructs.ConvertToNumpyArray( fp, arr )
        fplist.append( arr )
    return np.asarray( fplist )

res =calc_fp_arr(drugs)

kmeans = KMeans(init="random", n_clusters=nClusters, n_init=10, max_iter=400, random_state=42)
kmeans.fit(res)

min_distance = {}
cluster_center = {}
centeroid = {}
for i in range(nClusters):
    min_distance[i] = 999999
    cluster_center[i] = None
    centeroid[i] = kmeans.cluster_centers_[i]

for idx,fp in enumerate(res):
    cluster = kmeans.labels_[idx]
    print(fp.shape, centeroid[cluster].shape)
    distance = euclidean_distances([fp], [centeroid[cluster]])
    if distance < min_distance[cluster]:
        min_distance[cluster] = distance
        cluster_center[cluster] = idx


cluster_centers = cluster_center.values()

sd_writer = Chem.SDWriter(sys.argv[2])
for idx,mol in enumerate(drugs):
    if idx in cluster_centers:
        mol.SetProp("ClusterCenter","1")
    else:
        mol.SetProp("ClusterCenter","0")
    mol.SetProp("KMeanCluster", "%d"%kmeans.labels_[idx])
    sd_writer.write(mol)



