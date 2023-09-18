#!/usr/bin/env python

from openeye.oechem import *
from openeye.oemolprop import *
from openeye.oedepict import *
import sys, os,  glob, random,operator
import rdkit
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
from rdkit import DataStructs
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from datetime import datetime
import re
import time


if __name__ == "__main__":
    if len(sys.argv)!=4:
        print ("%s input.sdf numDiversityMols output.sdf"%sys.argv[0])
        sys.exit(1)
    molList = []
    ifs = oemolistream()
    ifs.open(sys.argv[1])
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        newmol = OEGraphMol(mol)
        molList.append(newmol)

    millis = int(round(time.time() * 1000))
    seed = int(str(millis)[-4:])

    ms = []
    dict = {}
    for oemol in molList:
        smiles = OEMolToSmiles(oemol)
        rdmol = Chem.MolFromSmiles(smiles)
        if rdmol is not None:
            mol_name = OEGetSDData(oemol,"IDNUMBER")
            rdmol.SetProp("_Name",mol_name)
            ms.append(rdmol)
            if mol_name in dict:
                print("Name collided!",file=sys.stderr)
            dict[mol_name] = oemol
        #    while ms.count(None): ms.remove(None)
    fps = [GetMorganFingerprint(x,3) for x in ms]
    def distij(i,j,fps=fps):
        dis = 1 - DataStructs.DiceSimilarity(fps[i], fps[j])
        return dis
    nfps = len(fps)
    picker = MaxMinPicker()
    pickIndices = picker.LazyPick(distij,nfps,int(sys.argv[2]),seed=seed)
    #print (list(pickIndices))
    picks = [dict[ms[x].GetProp("_Name")] for x in pickIndices]
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_SDF)
    for p in picks:
        OEWriteMolecule(ofs,p)
    ofs.close()

