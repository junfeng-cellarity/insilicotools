#!/usr/bin/env python

from openeye.oechem import *
import sys,os
# > <BIO-NUMBER>
# Unknown

if __name__=="__main__":
    if len(sys.argv)!=2:
        print "Usage:%s input.sdf"%sys.argv[0]
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        count_bio = 0
        count_missing = 0
        while OEReadMolecule(ifs,mol):
            bionumber = OEGetSDData(mol,"BIO-NUMBER")
            if bionumber is not None and len(bionumber) > 0 and bionumber!="Unknown":
                count_bio += 1
            else:
                count_missing += 1
        ifs.close()

        print "%d compounds have bionumber"%count_bio
        print "%d compounds do not have bionumber"%count_missing
