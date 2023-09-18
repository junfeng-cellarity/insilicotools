#!/usr/bin/env python
from openeye.oechem import *
import os, sys

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Usage:%s input.sdf input.csv output.sdf"
        print "The first line of input.csv should contain data tags"
    else:
        moldict = {}
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            mol_name = mol.GetTitle().decode("utf-8-sig").encode("utf-8")
            moldict[mol_name] = OEGraphMol(mol)
        ifs.close()

        csvfile = open(sys.argv[2],"r")
        lines = csvfile.read().splitlines()
        tags = lines[0].split(",")
        dict = {}
        for id,tag in enumerate(tags):
            dict[id] = tag

        ofs = oemolostream()
        ofs.open(sys.argv[3])
        for line in lines[1:]:
            args = line.split(",")
            mol = moldict[args[0]]
            OEClearSDData(mol)
            for id,arg in enumerate(args[1:]):
                tagName = dict[id+1]
                print tagName
                if len(arg.strip())==0:
                    continue
                else:
                    OESetSDData(mol,tagName,arg)
            OEWriteMolecule(ofs,mol)
        ofs.close()
        csvfile.close()
