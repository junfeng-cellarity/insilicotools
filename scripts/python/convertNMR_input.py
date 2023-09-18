#!/usr/bin/env python

import os, sys
from openeye.oechem import *
def get_atom_desc(atm):
    return "%s"%(atm.GetName())

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print ("Usage:%s input.pdb output.fmt"%sys.argv[0])
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        if OEReadMolecule(ifs,mol):
            residue = None
            print (mol.NumAtoms())
            for atom in mol.GetAtoms():
                residue = OEAtomGetResidue(atom)
                break
            n = 0
            for bond in mol.GetBonds():
                #bond = OEBondBase()
                if not bond.IsInRing():
                    bgn = bond.GetBgn()
                    end = bond.GetEnd()
                    #bgn = OEAtomBase()
                    #end = OEAtomBase()
                    if not bgn.IsHydrogen() and not end.IsHydrogen() and bgn.GetDegree()>1 and end.GetDegree()>1:
                        n = n+1
                        #print(bgn.GetName(),end.GetName())
            print("RESIDUE   %3s %5d %5d %4d %4d" % (residue.GetName(),n,mol.NumAtoms()+15,4,73))
            n = 1
            for bond in mol.GetBonds():
                # bond = OEBondBase()
                if not bond.IsInRing():
                    bgn = bond.GetBgn()
                    end = bond.GetEnd()
                    if not bgn.IsHydrogen() and not end.IsHydrogen():
                        atom1 = None
                        for atm in bgn.GetAtoms():
                            if atm.GetName() != end.GetName():
                                atom1 = atm
                                break
                        atom2 = None
                        for atm in end.GetAtoms():
                            if atm.GetName() != bgn.GetName():
                                atom2 = atm
                                break
                        if atom1 is not None and atom2 is not None:
                            name1 = atom1.GetName()
                            name2 = atom2.GetName()
                            print("%4d T%03d %5d%5d %9.4f  %-4s %-4s %-4s %-4s"%(n, n, 0, 0, 0.0, name1, bgn.GetName(), end.GetName(), name2))
                            n += 1
        ifs.close()
