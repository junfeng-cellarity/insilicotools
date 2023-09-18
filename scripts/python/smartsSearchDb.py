#!/usr/bin/env python
__author__ = 'jfeng1'
from openeye.oechem import *
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: %s pattern.sma database.sdf.gz"%sys.argv[0])
        sys.exit(1)
    smarts = open(sys.argv[1]).read().strip()
    p = OESubSearch()
    p.Init(smarts)

    ofs = oemolostream()
    ofs.open("hits.sdf")
    ofs.SetFormat(OEFormat_SDF)

    ifs = oemolistream()
    ifs.open(sys.argv[2])
    mol = OEGraphMol()
    n = 0
    while OEReadMolecule(ifs,mol):
        targetMol = OEGraphMol(mol)
        OEPrepareSearch(targetMol,p)
        if p.SingleMatch(targetMol):
            n += 1
            print("found "+str(n)+" hit(s).",file=sys.stderr)
            OEWriteMolecule(ofs,mol)
    ifs.close()
    ofs.close()
