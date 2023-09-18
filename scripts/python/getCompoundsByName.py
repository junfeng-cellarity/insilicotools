#!/usr/bin/env python
from openeye.oechem import *
from openeye.oeiupac import *
import os,sys
import re

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Usage:%s input.txt output.sdf"%sys.argv[0]
    else:
        nameFile = sys.argv[1]
        sdfFile = sys.argv[2]
        ofs = oemolostream()
        ofs.open(sdfFile)
        for line in open(nameFile,"r").read().splitlines():
            if len(line.strip())==0:
                continue
            if line.startswith("#"):
                continue
            args = line.split("\t")
            id = args[0]
            name = args[1]
#            name = re.sub(r'\(.*\)',"",args[1]).strip()
            if len(name) <= 3:
                print >> sys.stderr, "Failed to parse"
                print line
                continue
            else:
                mol = OEGraphMol()
                result = OEParseIUPACName(mol,name)
                if result:
                    OESetSDData(mol,"ID",id)
                    mol.SetTitle(name)
                    OEWriteMolecule(ofs,mol)
                else:
                    print >> sys.stderr, "Failed to parse"
                    print line
        ofs.close()