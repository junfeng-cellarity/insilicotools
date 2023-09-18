#!/usr/bin/env python
from openeye.oechem import *
import sys

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print ("Usage: %s pattern.sdf database.sdf core_hitnumber.sdf hit_structures.sdf"%sys.argv[0])
        sys.exit(1)

    pattern = []
    ifs = oemolistream()
    ifs.open(sys.argv[1])
    mol = OEGraphMol()
    idx = 0
    nameDict = {}
    molDict = {}
    while OEReadMolecule(ifs,mol):
        if mol.NumAtoms()>0:
            mol_name = mol.GetTitle()
            smarts = OEMolToSmiles(mol)
            p = OESubSearch(smarts)
            pattern.append(p)
            nameDict[idx] = mol_name
            molDict[idx] = OEGraphMol(mol)
            idx += 1
    ifs.close()

    matchResult = {}

    ifs = oemolistream()
    ifs.open(sys.argv[2])
    while OEReadMolecule(ifs,mol):
        for idx,p in enumerate(pattern):
            targetMol = OEGraphMol(mol)
            OEPrepareSearch(targetMol,p)
            if p.SingleMatch(targetMol):
                if idx not in matchResult:
                    matchResult[idx] = []
                matchResult[idx].append(targetMol)
    ifs.close()



    ofs = oemolostream()
    ofs.open(sys.argv[4])
    ofs.SetFormat(OEFormat_SDF)

    ofs2 = oemolostream()
    ofs2.open(sys.argv[3])
    ofs2.SetFormat(OEFormat_SDF)
    for idx,p in enumerate(pattern):
        mol = molDict[idx]
        numHits = 0
        if idx in matchResult:
            numHits = len(matchResult[idx])
        OESetSDData(mol,"HIT","%d"%numHits)
        OEWriteMolecule(ofs, mol)
        if idx in matchResult:
            for m in matchResult[idx]:
                OESetSDData(m,"HIT","%d"%idx)
                OESetSDData(m,"CoreName",nameDict[idx])
                OESetSDData(m,"CoreSmiles",OEMolToSmiles(molDict[idx]))
                OEWriteMolecule(ofs2,m)
    ofs2.close()
    ofs.close()
