
# coding: utf-8

# In[4]:


# import OE tools and data
#import oenotebook as oenb
from openeye.oechem import *
from openeye.oegraphsim import *

# In[5]:


import matplotlib.pyplot as plt
import numpy as np
#import mpld3
# get_ipython().magic(u'matplotlib inline')
import seaborn as sns
sns.set_context("talk")
colors = sns.color_palette()
# sns.palplot(colors)
do_tooltips = False

descriptor_dict = {
    "Circular":OEFPType_Circular,
    "Lingo":OEFPType_Lingo,
    "MACSS":OEFPType_MACCS166,
    "Path":OEFPType_Path,
    "Tree":OEFPType_Tree,
}

fpList = []
ifs = oemolistream()
ifs.open("/Users/jfeng1/BiogenDB/EVotec/evotec_descriptors.sdf")
mol = OEGraphMol()
molList = []
while OEReadMolecule(ifs,mol):
    fp_array = []
    fingerprint = OEFingerPrint()
    OEMakeFP(fingerprint,mol,descriptor_dict["Circular"])
    for i in range(0,fingerprint.GetSize()):
        if fingerprint.IsBitOn(i):
            fp_array.append(1)
        else:
            fp_array.append(0)
    fpList.append(fp_array)
    molList.append(OEGraphMol(mol))
ifs.close()



    
    


# In[6]:


# from sklearn.manifold import TSNE
# two_tsne = TSNE(metric="jaccard")
# Y = two_tsne.fit_transform(np.array(fpList))
from MulticoreTSNE import MulticoreTSNE as TSNE
import numpy as np
tsne = TSNE(n_jobs=8,metric="jaccard")
Y = tsne.fit_transform(np.array(fpList))
f = plt.figure()
sc = plt.scatter(Y[:,0], Y[:, 1], s=100, alpha=0.5,label="HTS")
plt.show()




