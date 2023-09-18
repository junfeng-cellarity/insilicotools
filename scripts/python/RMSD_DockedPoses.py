#!/usr/bin/env python
from openeye.oechem import *
import sys,os

def getMatchCount(match):
    n = 0
    for atm in match.GetAtoms():
        n += 1
    return n

if __name__ == "__main__":
    if len(sys.argv)!=4:
        print "Usage:%s template.sdf allPoses.sdf output.sdf"%sys.argv[0]
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        templateMol = OEGraphMol()
        OEReadMolecule(ifs,templateMol)
        ifs.close()

        mol = OEGraphMol()
        ifs = oemolistream()
        ifs.open(sys.argv[2])
        molList = []
        while OEReadMolecule(ifs,mol):
            molList.append(OEGraphMol(mol))
        ifs.close()

        mcs = OEMCSSearch(templateMol,OEExprOpts_DefaultAtoms,OEExprOpts_DefaultBonds,OEMCSType_Exhaustive)
        mcs.SetMCSFunc(OEMCSMaxAtoms())
        mcs.SetMinAtoms(10)
        unique = True

        ofs = oemolostream()
        ofs.open(sys.argv[3])
        for target in molList:
            all_matches = mcs.Match(target, unique)
            best_rmsd = 999 #any big number will do
            for count,match in enumerate(all_matches):
                rmsd = OERMSD(mcs.GetPattern(), target, match)
                if rmsd < best_rmsd:
                    best_rmsd = rmsd
            OESetSDData(target,"RMSD","%5.2f"%best_rmsd)
            OEWriteMolecule(ofs,target)
        ofs.close()
