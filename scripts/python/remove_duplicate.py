#!/usr/bin/env python
from openeye.oechem import *
from openeye.oeiupac import *
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print ("Remove duplicate molecules from input")
        print ("Usage:%s input.sdf output.sdf")
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        dict = {}
        n = 0
        ofs = oemolostream()
        ofs.open(sys.argv[2])
        while OEReadMolecule(ifs,mol):
            newmol=OEGraphMol(mol)
            OEDeleteEverythingExceptTheFirstLargestComponent(newmol)
            smiles = OEMolToSmiles(newmol)
            if smiles in dict:
                n += 1
                print ("%d duplicates found."%n)
            else:
                dict[smiles] = 1
                # name = OECreateIUPACName(mol)
                # if "BLAH" not in name:
                #     OESetSDData(mol,"IUPAC_NAME",name)
                # else:
                #     OESetSDData(mol,"IUPAC_NAME", "")
                OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()
