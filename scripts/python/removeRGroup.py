#!/usr/bin/env python
from openeye.oechem import *
from openeye.oeiupac import *
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "Remove Rgroup from input"
        print "Usage:%s input.sdf output.sdf"
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        ofs = oemolostream()
        ofs.open(sys.argv[2])
        ofs.SetFormat(OEFormat_SDF)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            mol1 = OEGraphMol(mol)
            atomToBeRemoved = []
            bondToBeRemoved = []
            for atom in mol1.GetAtoms():
                if atom.GetAtomicNum()==0:
                    atom.SetAtomicNum(1)
                    if atom.GetExplicitDegree()==1:
                        for nbr in atom.GetAtoms():
                            if nbr.GetExplicitDegree()==1:
                                atomToBeRemoved.append(nbr)
                                atomToBeRemoved.append(atom)
                                bondToBeRemoved.append(mol1.GetBond(atom,nbr))
            for atom in atomToBeRemoved:
                mol1.DeleteAtom(atom)
            for bond in bondToBeRemoved:
                mol1.DeleteBond(bond)
            OESuppressHydrogens(mol1)
            molName = OEGetSDData(mol1,"Scaffold_ID")
            mol1.SetTitle(molName)
            OEWriteMolecule(ofs,mol1)

        ifs.close()
        ofs.close()