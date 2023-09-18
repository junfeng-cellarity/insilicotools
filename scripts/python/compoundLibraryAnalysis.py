#!/usr/bin/env python
from openeye.oechem import *
from openeye.oegraphsim import *
import sys, os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

if len(sys.argv)!=4:
    print ("Usage:%s cmpdLibrary1.sdf cmpdLibrary2.sdf output.sdf"%sys.argv[0])
else:
    ifs1 = oemolistream()
    ifs1.open(sys.argv[1])
    molList = []
    fplist1 = []
    mol = OEGraphMol()
    while OEReadMolecule(ifs1,mol):
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, OEFPType_Circular)
        fplist1.append(fp)
        molList.append(OEGraphMol(mol))
    ifs1.close()

    ifs2 = oemolistream()
    ifs2.open(sys.argv[2])
    fplist2 = []
    neighborMolList = []
    while OEReadMolecule(ifs2,mol):
        fp = OEFingerPrint()
        OEMakeFP(fp, mol, OEFPType_Circular)
        fplist2.append(fp)
        neighborMolList.append(OEGraphMol(mol))
    ifs2.close()

    ofs = oemolostream()
    ofs.open(sys.argv[3])


    dist = []
    for id,fp1 in enumerate(fplist1):
        mol = molList[id]
        best_cmpd_id = None
        best_tanimoto = 0
        for id2,fp2 in enumerate(fplist2):
            tanimoto = OETanimoto(fp1,fp2)
            if tanimoto > best_tanimoto:
                best_cmpd_id = id2
                best_tanimoto = tanimoto
        neighbor = neighborMolList[best_cmpd_id]
        smiles = "%s %s"%(OEMolToSmiles(neighbor),neighbor.GetTitle())
        OESetSDData(mol,"Best Tanimoto","%5.2f"%best_tanimoto)
        OESetSDData(mol,"Best Neighbor",smiles)
        OEWriteMolecule(ofs,mol)
        dist.append(best_tanimoto)
    ofs.close()
    # stddev = np.std(np.array(dist))
    # mean = np.mean(np.array(dist))
    # n, bins, patches = plt.hist(np.array(dist),20, normed=True, facecolor='green', alpha=0.5)
    # y = mlab.normpdf(bins, mean, stddev)
    # plt.plot(bins, y, 'r--')
    # plt.xlabel('Tanimoto')
    # plt.ylabel('Probability')
    # plt.title('Histogram of Similarity Distribution:')
    #
    # # Tweak spacing to prevent clipping of ylabel
    # # plt.subplots_adjust(left=0.15)
    # plt.show()

