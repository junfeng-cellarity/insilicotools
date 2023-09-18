#!/usr/bin/env python
from openeye.oechem import *
import sys

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print "Usage: %s input.sdf output.sdf tag isBiggerBetter"%sys.argv[0]
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])

        ofs = oemolostream()
        ofs.open(sys.argv[2])

        tag  = sys.argv[3].strip()
        isBiggerBetter = int(sys.argv[4])>0

        dict = {}
        mol = OEGraphMol()
        name_list = []
        while OEReadMolecule(ifs,mol):
            name = mol.GetTitle().split("_")[0]
            if name not in name_list:
                name_list.append(name)
            if not dict.has_key(name):
                dict[name] = OEGraphMol(mol)
            else:
                value1 = float(OEGetSDData(mol,tag))
                value2 = float(OEGetSDData(dict[name],tag))
                if isBiggerBetter and value1>value2:
                    dict[name] = OEGraphMol(mol)
                if not isBiggerBetter and value1<value2:
                    dict[name] = OEGraphMol(mol)
        for name in name_list:
            mol = dict[name]
            mol.SetTitle(name)
            OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()
