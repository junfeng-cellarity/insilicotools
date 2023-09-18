#!/usr/bin/env python
import sys
from openeye.oechem import *

if __name__=="__main__":
    num = int(sys.argv[2])
    ifs = oemolistream()
    ifs.open(sys.argv[1])
    mol = OEGraphMol()
    OEReadMolecule(ifs,mol)
    ifs.close()
    # OEPerceiveResidues(mol)
    resList = []
    for atm in mol.GetAtoms():
        res = OEAtomGetResidue(atm)
        print res.GetResidueNumber()
        # res.SetResidueNumber(res.GetResidueNumber()+num)

    ofs = oemolostream()
    ofs.SetFormat(OEFormat_PDB)
    ofs.openstring()
    OEWriteMolecule(ofs,mol)
    print ofs.GetString()
    ofs.close()

