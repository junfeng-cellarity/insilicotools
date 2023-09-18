#!/usr/bin/env python
from openeye.oechem import *
import os
import sys

if __name__ == "__main__":
    if len(sys.argv)!=2:
        print "Usage:%s input.sdf"%(sys.argv[0])
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])

        mol = OEGraphMol()
        dict = []
        while OEReadMolecule(ifs,mol):
            OEPerceiveChiral(mol)
            OEMDLPerceiveBondStereo(mol)
            molName = OEGetSDData(mol,"Alias_ID")
            for bond in mol.GetBonds():
                # bond = OEBondBase()
                if bond.HasStereoSpecified(OEBondStereo_Wedge):
                    print molName
                    if bond.GetBgn().IsNitrogen() or bond.GetEnd().IsNitrogen:
                        if molName not in dict:
                            dict.append(molName)
                            print molName

        ifs.close()
