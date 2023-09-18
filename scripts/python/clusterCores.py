#!/usr/bin/env python
from openeye.oechem import *
from openeye.oemedchem import *
import sys, os

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print "Split sdf file into different bins based on cores, each structure should only contain one core."
        print "%s input.sdf output.sdf"%sys.argv[0]
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        ofs = oemolostream()
        ofs.open(sys.argv[2])
        dict = {}
        coreKeys = []
        while OEReadMolecule(ifs,mol):
            newmol = OEGraphMol(mol)
            oechem.OETheFunctionFormerlyKnownAsStripSalts(mol)
            frameworks = OEGetBemisMurcko(mol, OERegionType_Framework)
            for framework in frameworks:
                frag = OEGraphMol()
                oechem.OESubsetMol(frag,newmol,framework)
                key = OEMolToSmiles(frag)
                if not dict.has_key(key):
                    dict[key] = []
                    dict[key].append(newmol)
                    coreKeys.append(key)
                else:
                    dict[key].append(newmol)
                break
        ifs.close()

        sorted(coreKeys)

        for core in coreKeys:
            for mol in dict[core]:
                OESetSDData(mol,"coreSmiles",core)
                OESetSDData(mol,"coreId","%d"%coreKeys.index(core))
                OEWriteMolecule(ofs,mol)
        ofs.close()

        ofs.close()
