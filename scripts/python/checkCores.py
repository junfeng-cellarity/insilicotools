#!/usr/bin/env python

import glob
import os
from openeye.oechem import *

sdf = "/Users/jfeng1/WatchDox/WuXi/WuXi_53cores.sdf"
directorylist = glob.glob("/Users/jfeng1/WatchDox/WuXi/WuXi_53_Cores_*")
dict = []
for d in directorylist:
    files = glob.glob("%s/*.sdf"%d)
    for f in files:
        fname = os.path.basename(f)
        args = fname.split("_", 1)
        if args[0]=="BIO":
            bio_number = "BIO-%s"%(args[1].split("_")[0])
        else:
            bio_number = args[0]
        if bio_number not in dict:
            dict.append(bio_number)
print len(dict)

ifs = oemolistream()
ifs.open(sdf)
mol = OEGraphMol()
while OEReadMolecule(ifs,mol):
    bio = mol.GetTitle()
    if bio not in dict:
        print bio
ifs.close()