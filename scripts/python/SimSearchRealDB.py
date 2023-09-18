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
import glob
from scipy.stats import norm
from heapq import heapify, heappush, heappushpop
class MaxHeap():
    def __init__(self, top_n):
        self.h = []
        self.length = top_n
        heapify( self.h)

    def add(self, element):
        if len(self.h) < self.length:
            heappush(self.h, element)
        else:
            heappushpop(self.h, element)

    def getTop(self):
        return sorted(self.h, reverse=True)

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def get_best_hits(smiles_file, probe_fp, hit_list, max_hits):
    progressbar = ProgressBar(maxval=file_len(smiles_file)).start()
#    progressbar = ProgressBar(maxval=200000).start()
    file = open(smiles_file, "r")
    result_heap = MaxHeap(max_hits)
    id = 0
    for line in file:
        smiles,name = line.strip().split()
        mol = Chem.MolFromSmiles(smiles)
#       fp = Chem.AllChem.GetMorganFingerprint(mol,3)
        fp = Chem.AllChem.GetAtomPairFingerprint(mol)
        similarity = DataStructs.TanimotoSimilarity(probe_fp,fp)
        result_heap.add((similarity,smiles,name))
        progressbar.update(id)
        id = id+1
    file.close()
    hit_list.append(result_heap.getTop())
    progressbar.finish()

if len(sys.argv)!=3:
    print ("Usage:%s probe.sdf output.sdf"%sys.argv[0])
else:

    fplist1 = []
    file1 = Chem.SDMolSupplier(sys.argv[1])
    sameLibrary = False
    if sys.argv[1] == sys.argv[2]:
        sameLibrary = True
    probe_fp = None
    for mol in file1:
        try:
            probe_fp = AllChem.GetAtomPairFingerprint(mol)
            break
        except:
            continue
    max_hits = 500
#    files = glob.glob("/home/jfeng/Database/RealDB_2.0/MolTrunk*.smi")
    files = glob.glob("/home/jfeng/Database/RealDB_2.0/realdb_part*.smi")
    manager = multiprocessing.Manager()
    hit_list = manager.list()
    jobs = []
    for idx,file in enumerate(files):
        p = multiprocessing.Process(target=get_best_hits, args=(file,probe_fp,hit_list,max_hits))
        jobs.append(p)
        p.start()
    for p in jobs:
        p.join()

    print(hit_list)
    sd_writer = Chem.SDWriter(sys.argv[2])
    for mylist in hit_list:
        for sim,smiles,name in mylist:
            mol = Chem.MolFromSmiles(smiles)
            mol.SetProp("_Name",name)
            mol.SetProp("Tanimoto","%5.2f"%sim)
            sd_writer.write(mol)
    sd_writer.close()

