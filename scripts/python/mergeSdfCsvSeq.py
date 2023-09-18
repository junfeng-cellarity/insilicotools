#!/usr/bin/env python
from openeye.oechem import *
import sys
if __name__=="__main__":
    if len(sys.argv)!=4:
        print "Usage:%s input.sdf input.csv output.sdf"%sys.argv[0]
        print "merge sdf with input.csv in order, first line of csv file contains tag names"
    else:
        csvDict = {}
        csv = open(sys.argv[2],"r")

        colNames = []
        for lineno,line in enumerate(csv.read().splitlines()):
            if lineno == 0:
                colNames = line.split(",")
            else:
                args = line.split(",")
                csvDict[lineno] = {}
                for id,arg in enumerate(args):
                    csvDict[lineno][colNames[id].replace("\"","")] = args[id].replace("\"","")
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        oemol = OEGraphMol()
        ofs = oemolostream()
        ofs.open(sys.argv[3])
        molid = 1
        while OEReadMolecule(ifs,oemol):
            mol = OEGraphMol(oemol)
            tagDict = csvDict[molid]
            for key in tagDict.keys():
                OESetSDData(mol,key,tagDict[key])
            OEWriteMolecule(ofs,mol)
            molid += 1
        ofs.close()
        ifs.close()

