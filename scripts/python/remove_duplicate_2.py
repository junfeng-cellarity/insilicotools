#!/usr/bin/env python
__author__ = 'jfeng1'
from openeye.oechem import *
from openeye.oeiupac import *
import sys

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print ("Remove duplicate molecules contained in second file from first file")
        print ("Usage:%s input.sdf exclude.sdf output.sdf duplicate.sdf")
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[2])
        mol = OEGraphMol()
        dict = {}
        while OEReadMolecule(ifs,mol):
            OEDeleteEverythingExceptTheFirstLargestComponent(mol)
            dict[OEMolToSmiles(mol)] = 1
        ifs.close()

        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        n = 0
        ofs = oemolostream()
        ofs.open(sys.argv[3])
        ofs.SetFormat(OEFormat_SDF)

        ofs2 = oemolostream()
        ofs2.open(sys.argv[4])
        ofs2.SetFormat(OEFormat_SDF)
        while OEReadMolecule(ifs,mol):
            newmol = OEGraphMol(mol)
            OEDeleteEverythingExceptTheFirstLargestComponent(newmol)
            smiles = OEMolToSmiles(newmol)
            if smiles in dict:
                n += 1
                print (sys.stderr, "%d duplicates found."%n)
                OEWriteMolecule(ofs2,mol)
            else:
                OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()
        ofs2.close()
