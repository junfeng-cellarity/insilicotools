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

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i+n]

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def get_best_hits(smiles_file, pattern, hit_list):
    progressbar = ProgressBar(maxval=file_len(smiles_file)).start()
#    progressbar = ProgressBar(maxval=200000).start()
    file = open(smiles_file, "r")
    id = 0
    for line in file:
        smiles,name = line.strip().split()
        mol = Chem.MolFromSmiles(smiles)
        if mol.HasSubstructMatch(pattern):
            hit_list.append((mol,name))
        progressbar.update(id)
        id = id+1
    file.close()
    progressbar.finish()

if len(sys.argv)!=3:
    print ("Usage:%s probe.smarts output.sdf"%sys.argv[0])
else:
    smarts = open(sys.argv[1],"r").read().strip()
    pattern = Chem.MolFromSmarts(smarts)
    #files = glob.glob("/home/jfeng/Database/RealDB_2.0/test/MolTrunk*.smi")
    files = glob.glob("/home/jfeng/Database/RealDB_2.0/realdb_part*.smi")
    manager = multiprocessing.Manager()
    hit_list = manager.list()
    jobs = []
    for idx,file in enumerate(files):
        p = multiprocessing.Process(target=get_best_hits, args=(file,pattern,hit_list))
        jobs.append(p)
        p.start()
    for p in jobs:
        p.join()

    print("%d hits found."%len(hit_list))
    sd_writer = Chem.SDWriter(sys.argv[2])
    for mol,name in hit_list:
        mol.SetProp("_Name",name)
        sd_writer.write(mol)
    sd_writer.close()

