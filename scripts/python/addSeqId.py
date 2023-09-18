#!/usr/bin/env python
from openeye.oechem import *
import sys
import os

if __name__=="__main__":
    if len(sys.argv)!=3:
        print "Add a SeqId tag to the sdf file"
        print "Usage:%s input.sdf output.sdf"%(os.path.basename(sys.argv[0]))
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])

        ofs = oemolostream()
        ofs.open(sys.argv[2])

        mol = OEGraphMol()
        id = 0
        while OEReadMolecule(ifs, mol):
            OESetSDData(mol,"SeqId","%d"%id)
            OEWriteMolecule(ofs,mol)
            id += 1
        ifs.close()
        ofs.close()
