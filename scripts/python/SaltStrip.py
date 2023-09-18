#!/usr/bin/env python
__author__ = 'jfeng1'
import sys
from openeye.oechem import *

ifs = oemolistream()
ifs.open(sys.argv[1])
mol = OEGraphMol()

ofs = oemolostream()
ofs.open(sys.argv[2])

while OEReadMolecule(ifs,mol):
    numparts, partlist = OEDetermineComponents(mol)
    if numparts>1:
        pred = OEPartPredAtom(partlist)
        maxAtoms = 0
        partId = -1
        for i in range(1, numparts + 1):
            pred.SelectPart(i)
            partmol = OEGraphMol()
            OESubsetMol(partmol, mol, pred)
            if partmol.NumAtoms() > maxAtoms:
                maxAtoms = partmol.NumAtoms()
                partId = i

        smiles = []
        for i in range(1,numparts+1):
            if i == partId:
                continue
            pred.SelectPart(i)
            partmol = OEGraphMol()
            OESubsetMol(partmol, mol, pred)
            smiles.append(OEMolToSmiles(partmol))

        newmol = OEGraphMol()
        pred.SelectPart(partId)
        OESubsetMol(newmol,mol,pred)
        OESetSDData(newmol,"Smiles of Removed Salts",".".join(smiles))
        OEWriteMolecule(ofs,newmol)
    else:
        OEWriteMolecule(ofs,mol)
ifs.close()
ofs.close()
