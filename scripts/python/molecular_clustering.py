#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
from rdkit import Chem

if len(sys.argv)!=4:
    print ("Usage: %s input.sdf cutoff output.sdf")
    sys.exit(1)


ifs = oemolistream()
ifs.open(sys.argv[1])
mol = OEGraphMol()
simMatrix = {}
mols = []
while OEReadMolecule(ifs,mol):
    mols.append(OEGraphMol(mol))
ifs.close()


# In[2]:


fingerprint_db = OEFPDatabase(OEFPType_Circular)

emptyfp = OEFingerPrint()
emptyfp.SetFPTypeBase(fingerprint_db.GetFPTypeBase())

for idx,mol in enumerate(mols):
    fingerprint_db.AddFP(mol)


# In[ ]:


cutoff = float(sys.argv[2])
def getNeighbors(idx):
    neighbors = []
    mol = mols[idx]
    numFps = fingerprint_db.NumFingerPrints()
    fingerprint_db.SetCutoff(cutoff)
    scores = fingerprint_db.GetScores(mol,0,numFps)
    for score in scores:
        if score.GetIdx()!=idx:
            neighbors.append(score.GetIdx())
    return (idx,neighbors)

import time
start = time.time()
import progressbar
import multiprocessing
clusterDict = {}
marker_dict = {}
nProcessor = 5
pool = multiprocessing.Pool(nProcessor)
argList = []
for idx,mol in enumerate(mols):
    marker_dict[idx] = 0
    argList.append(idx)

progressbar = progressbar.ProgressBar(maxval=len(argList))
rs = pool.map_async(getNeighbors,argList,chunksize=1)
progressbar.start()
while True:
    if rs.ready():
        break
    else:
        progress = len(argList)-rs._number_left
        progressbar.update(progress)
        time.sleep(0.5)
result = rs.get()
for idx,nbrList in result:
    clusterDict[idx] = nbrList
sortedList = sorted(clusterDict,key=lambda k:len(clusterDict[k]),reverse=True)
#print sortedList
end = time.time()
print 
print (end-start)

clusters = []
cluster_0 = clusterDict[sortedList[0]]
cluster_0.append(sortedList[0])
clusters.append(cluster_0)
for x in cluster_0:
    marker_dict[x] = 1
    
for x in sortedList[1:]:
    if marker_dict[x] == 1:
        continue
    else:
        cluster = []
        cluster.append(x)
        marker_dict[x] = 1
        for a in clusterDict[x]:
            if marker_dict[a] ==0:
                cluster.append(a)
                marker_dict[a] = 1
        clusters.append(cluster)
print(len(clusters))
file = open("/Users/jfeng1/cluster_center.txt","w")
for idx,c in enumerate(clusters):
    print(OEMolToSmiles(mols[c[0]]),mols[c[0]].GetTitle(),file=file)
file.close()


ofs = oemolostream()
ofs.open(sys.argv[3])
for idx,c in enumerate(clusters):
    for x in c:
        OESetSDData(mols[x],"cluster_no","%d"%idx)
        OESetSDData(mols[x],"cluster_center","0")
    OESetSDData(mols[c[0]],"cluster_center","1")
    OEWriteMolecule(ofs,mols[c[0]])
ofs.close()

