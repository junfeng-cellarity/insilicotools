#!/usr/bin/env python
from openeye.oechem import *
import sys

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print "Usage: %s reaction.smirk reactant1.smi reactant2.smi ..."%sys.argv[0]
    else:
        smirkFile = open(sys.argv[1],"r")
        smirk = smirkFile.read().split("\n")[0].strip()
        print smirk
        # smirk = "[NH2:1]>>[Nh3+:1]"
        libraryGen = OELibraryGen()
        libraryGen.SetExplicitHydrogens(False)
        libraryGen.SetValenceCorrection(True)
        libraryGen.Init(smirk)
        n = libraryGen.NumReactants()
        if len(sys.argv) < n+2:
            print "Usage: %s reaction.smirk reactant1.smi reactant2.smi ..."%sys.argv[0]
        else:
            for i in range(n):
                ifs = oemolistream()
                ifs.open(sys.argv[i+2])
                reactant = OEGraphMol()
                while OEReadMolecule(ifs,reactant):
                    libraryGen.AddStartingMaterial(OEGraphMol(reactant),i)
            for product in libraryGen.GetProducts():
                print OEMolToSmiles(product), product.GetTitle()
