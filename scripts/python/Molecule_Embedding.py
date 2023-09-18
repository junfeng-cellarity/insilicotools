#!/usr/bin/env python
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_context("talk")

from openeye.oechem import *
from openeye.oegraphsim import *
import numpy as np
import scipy as sp
from scipy.spatial.distance import squareform
import progressbar
import bhtsne

class CachingMolDistance():

    def __init__(self):
        self.fpIdx = []
        self.nMols = 0
        self.DistArray = np.empty(0)

    def __repr__(self):
        return("Caching Distance: Pairwaise distances of {} molecules cached".format(self.nMols))

    def dist(self, fpA, fpB):
        return OETanimoto(fpA,fpB)

    def load(self, fp_iterable):
        size = len(fp_iterable)
        pbar = progressbar.ProgressBar(size).start()
        dv = []
        for i in range(size):
            pbar.update(i)
            for j in range(i+1,size):
                dv.append(self.dist(fp_iterable[i],fp_iterable[j]))
        self.DistArray = squareform(np.asarray(dv))


#molfile = "/Users/jfeng1/Datasets/Evotec/biogen_in_stock_april_11_2018.sdf"
import sys
selected_file = sys.argv[2]
selected_dict = {}
selected_list = []

mol = OEGraphMol()
ifs = oemolistream()
ifs.open(selected_file)
while OEReadMolecule(ifs,mol):
    selected_dict[OEMolToSmiles(mol)] =1


molfile = sys.argv[1]
ifs = oemolistream()
ifs.open(molfile)

mol = OEGraphMol()
molList = []
fpList = []
count = 0

while OEReadMolecule(ifs,mol):
    smiles = OEMolToSmiles(mol)
    if smiles in selected_dict:
        selected_list.append(count)
    molList.append(OEGraphMol(mol))
    fp = OEFingerPrint()
    OEMakeFP(fp,mol,OEFPType_Circular)
    fpList.append(fp)
    count += 1
ifs.close()

selected_list = np.array(selected_list)

mol_dists = CachingMolDistance()
mol_dists.load(fpList)

# from sklearn.manifold import MDS
# mds = MDS(dissimilarity="precomputed", max_iter=1000)
# Y = mds.fit_transform(mol_dists.DistArray)
#
# # In[12]:
#
# f = plt.figure()
# sc = plt.scatter(Y[:, 0], Y[:, 1], s=100, alpha=0.5,label="All")
# plt.scatter(Y[selected_list,0],Y[selected_list,1], color = "red", alpha=0.1)
#
# plt.title("MDS")
# plt.text(0.8*max(Y[:,0]),0.95*max(Y[:,1]),"All",color="blue",size=24)
# plt.show()

#
# # ### Non-metric MDS
# #
# # Multidimensional scaling can also be performed without the assumption that this distances are a metric which satisfy the triangle inequality. The embedding attempts to preserve the order of the distances.
#
# # In[13]:
#
#
# nmds = MDS(dissimilarity="precomputed", max_iter=1000, metric=False)
# Y = nmds.fit_transform(mol_dists.DistArray)
# f = plt.figure()
# sc = plt.scatter(Y[:, 0], Y[:, 1], s=100, alpha=0.5,label="All")
#
# plt.title("nMDS")
# plt.text(0.8*max(Y[:,0]),0.95*max(Y[:,1]),"All",color="blue",size=24)
# plt.show()

# ## t-distributed Stochastic Neighbor Embedding
# 
# The **t-SNE** algorithm treats distances as joint gaussian distributions and minimizes the [KL divergence](https://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence) between the probabilities in the low-dimensional embedding and the original high-dimensional space.

# In[15]:


from sklearn.manifold import TSNE


# In[16]:

tsne = TSNE(metric="precomputed", learning_rate=150)
Y = tsne.fit_transform(mol_dists.DistArray)


# In[17]:


f = plt.figure()
plt.scatter(Y[:, 0], Y[:, 1], s=100, alpha=1.0,label="All")
plt.scatter(Y[selected_list,0],Y[selected_list,1], color = "red", alpha=0.3)

plt.title("t-SNE")
plt.text(0.8*max(Y[:,0]),0.95*max(Y[:,1]),"All",color="blue",size=24)
plt.show()

sys.exit(0)



# # ## Spectral Embedding
# #
# # Spectral embedding finds a low-dimensional representation through spectral decomposition on the laplacian of the affinity graph.
# from sklearn.manifold import SpectralEmbedding
# se = SpectralEmbedding(affinity="precomputed")
# Y = se.fit_transform(1.0 - mol_dists.DistArray)
# f = plt.figure()
# sc = plt.scatter(Y[:, 0], Y[:, 1], s=100, alpha=0.5,label="All")
# plt.text(0.8*max(Y[:,0]),0.95*max(Y[:,1]),"All",color="blue",size=24)
#
# plt.title("Spectral Embedding Two Targets")
# plt.show()