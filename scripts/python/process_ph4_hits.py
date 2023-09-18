#!/usr/bin/env python
from openeye.oechem import *
import sys
if len(sys.argv)!=3:
    print("%s input.sdf output.sdf"%sys.argv[0])
    sys.exit(1)
molecules = {}
mol_names = []
ifs = oemolistream()
ifs.open(sys.argv[1])
mol = OEGraphMol()
while OEReadMolecule(ifs,mol):
    mol_name = OEGetSDData(mol,"Molecule Name")
    if mol_name not in mol_names:
        mol_names.append(mol_name)
        mol.SetTitle(mol_name)
        molecules[mol_name] = OEGraphMol(mol)
    else:
        new_score = float(OEGetSDData(mol,"r_phase_PhaseScreenScore"))
        old_score = float(OEGetSDData(molecules[mol_name],"r_phase_PhaseScreenScore"))
        if new_score>old_score:
            molecules[mol_name] = OEGraphMol(mol)
ifs.close()

ofs = oemolostream()
ofs.open(sys.argv[2])
for mol_name in mol_names:
    mol = molecules[mol_name]
    OEWriteMolecule(ofs,mol)
ofs.close()
