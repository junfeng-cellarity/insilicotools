#!/usr/bin/env python

from openeye.oechem import *

csv = "/Users/jfeng1/Datasets/SoftFocus/OGA.csv"

dict = {}
lines = open(csv,"r").read().splitlines()
for lineno,line in enumerate(lines):
    if lineno == 0:
        continue
    args = line.split(",")
    dict[args[0]] = 1


sdf = "/Users/jfeng1/Datasets/SoftFocus/output.sdf"

ifs = oemolistream()
ifs.open(sdf)
mol = OEGraphMol()

ofs = oemolostream()
ofs.open("oga_output.sdf")
while OEReadMolecule(ifs,mol):
    abase_id = OEGetSDData(mol,"Abase_ID")
    if abase_id in dict:
        OEWriteMolecule(ofs,mol)
ofs.close()