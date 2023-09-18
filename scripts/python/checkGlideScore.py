#!/usr/bin/env python
import os,glob
from openeye.oechem import *

directory = "/Users/jfeng1/TTBK1"
all_sdf = "ttbk1_result.sdf"

file = os.path.join(directory,all_sdf)

sdfs = glob.glob(os.path.join(directory,"*selected*.sdf"))

dict = {}

for sdf in sdfs:
    mol = OEGraphMol()
    ifs = oemolistream()
    ifs.open(sdf)
    while OEReadMolecule(ifs,mol):
        if OEHasSDData(mol, "r_i_glide_gscore"):
            glide_score = float(OEGetSDData(mol,"r_i_glide_gscore"))
            if dict.has_key(mol.GetTitle()):
                if glide_score < dict[mol.GetTitle()]:
                    dict[mol.GetTitle()] = glide_score
            else:
                dict[mol.GetTitle()] = glide_score
    ifs.close()

ofs = oemolostream()
ofs.open(os.path.join(directory,"ttbk1_result_glidescore.sdf"))
ifs = oemolistream()
ifs.open(file)
mol = OEGraphMol()
while OEReadMolecule(ifs,mol):
    glide_score = 99999
    if dict.has_key(mol.GetTitle()):
        OESetSDData(mol,"has_glide_score","YES")
        OESetSDData(mol,"glide_score","%5.3f"%dict[mol.GetTitle()])
    else:
        OESetSDData(mol,"has_glide_score","NO")
    OEWriteMolecule(ofs,mol)
ifs.close()
ofs.close()


