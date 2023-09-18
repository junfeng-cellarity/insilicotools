#!/usr/bin/env python
import os,sys
import schrodinger
from schrodinger.protein.getpdb import get_pdb

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: %s pdb_list.txt"%(sys.argv[0]))
    else:
        for line in open(sys.argv[1]):
            pdb_name = line.strip()
            get_pdb(pdb_name)
