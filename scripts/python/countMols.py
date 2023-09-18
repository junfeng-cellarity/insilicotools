#!/usr/bin/env python

from openeye.oechem import *
import sys
if __name__ == "__main__":
    if len(sys.argv)!=2:
        print ("Usage:%s input.sdf"%sys.argv[0])
    molDb = OEMolDatabase(sys.argv[1])
    print(molDb.GetMaxMolIdx())

