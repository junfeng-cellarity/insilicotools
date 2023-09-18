#!/usr/bin/env python
from openeye.oechem import *
import os
import sys

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print "Usage:%s input.sdf output.sdf"%(sys.argv[0])
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])

        ofs = oemolostream()
        ofs.open(sys.argv[2])

        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            comment = OEGetSDData(mol, "Stereochemistry")
            print mol.GetTitle(),OEMDLGetParity(mol),comment
            if len(comment.strip()) == 0:
                OEMDLSetParity(mol,False)
            elif comment =="Mixture of diastereomers":
                OEMDLSetParity(mol,False)
            else:
                OEMDLSetParity(mol,True)
            OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()
