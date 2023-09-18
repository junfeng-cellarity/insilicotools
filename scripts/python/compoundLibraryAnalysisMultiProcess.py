#!/usr/bin/env python
import sys, os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from rdkit import Chem
import multiprocessing
from progressbar import *
from rdkit.Chem import AllChem
from rdkit import DataStructs
from scipy.stats import norm

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def get_best_and_mean_tanimoto(probes, targets, all_sim, max_dict, mean_dict, sameLibrary):
    progress = ProgressBar(maxval=len(probes)).start()
    for id, values in enumerate(probes):
        probe, id1 = values
        best_tanimoto = 0
        mean_tanimoto = 0
        n = 0
        #progress2 = ProgressBar(maxval=len(targets)).start()
        for id2,target in enumerate(targets):
            if sameLibrary and id1==id2:
                continue
            #progress2.update(id2)
            tanimoto = DataStructs.TanimotoSimilarity(probe,target)
            all_sim.append(tanimoto)
            mean_tanimoto += tanimoto
            n += 1
            if tanimoto > best_tanimoto:
                best_tanimoto = tanimoto
        mean_tanimoto = mean_tanimoto/n
        #progress2.finish()
        max_dict[id1]=best_tanimoto
        mean_dict[id1] = mean_tanimoto
        progress.update(id)
    progress.finish()

if len(sys.argv)!=4:
    print ("Usage:%s cmpdLibrary1.sdf cmpdLibrary2.sdf output.sdf"%sys.argv[0])
else:

    fplist1 = []
    file1 = Chem.SDMolSupplier(sys.argv[1])
    file2 = Chem.SDMolSupplier(sys.argv[2])
    sameLibrary = False
    if sys.argv[1] == sys.argv[2]:
        sameLibrary = True
    id = 0
    molDict = {}
    for mol in file1:
        try:
            fp = AllChem.GetMorganFingerprint(mol,3)
            fplist1.append((fp,id))
            molDict[id] = mol
            id = id + 1
        except:
            continue

    fplist2 = []
    for mol in file2:
        try:
            fp = AllChem.GetMorganFingerprint(mol,3)
            fplist2.append(fp)
        except:
            continue

    n_processors = 40
    trunk_size = int(len(fplist1)/n_processors)
    trunks = list(chunks(fplist1,trunk_size))
    manager = multiprocessing.Manager()
    max_dict = manager.dict()
    mean_dict = manager.dict()
    all_sim = manager.list()
    jobs = []
    for idx,trunk in enumerate(trunks):
        p = multiprocessing.Process(target=get_best_and_mean_tanimoto, args=(trunk,fplist2,all_sim, max_dict,mean_dict,sameLibrary))
        jobs.append(p)
        p.start()
    for p in jobs:
        p.join()


    sd_writer = Chem.SDWriter(sys.argv[3])
    ids = sorted(max_dict.keys())
    for id in ids:
        mol = molDict[id]
        mol.SetProp("MaxSim", str(max_dict[id]))
        mol.SetProp("MeanSim",str(mean_dict[id]))
        sd_writer.write(mol)
    sd_writer.close()
    narray = np.array(all_sim)
    print("Median:%f"%np.median(narray))
    print("N:%d"%len(all_sim))
    print(np.quantile(narray,0.95))


    # final_result = []
    # for dis in result_dict.values():
    #     final_result.append(dis)
    # sample = np.array(final_result)
    # stddev = np.std(sample)
    # mean = np.mean(sample)
    # weights = np.ones_like(sample) / (len(sample))
    # # x = np.linspace(0, 1, 100)
    # # plt.plot(x, norm.pdf(x, mean, stddev))
    # n, bins, patches = plt.hist(sample, 50, weights=weights, facecolor='green', alpha=0.5)
    # # y = norm.pdf(bins, mean, stddev)
    # # plt.plot(bins, y, 'r--')
    # plt.xlabel('Tanimoto')
    # plt.ylabel('Probability')
    # plt.title('Histogram of Similarity Distribution:')
    # plt.savefig("output.png")

