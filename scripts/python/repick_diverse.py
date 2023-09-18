#!/usr/bin/env python
import rdkit
from rdkit import Chem,DataStructs
import os,sys
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint
import numpy as np

if __name__ == "__main__":
    if len(sys.argv)!=5:
        print("Usage: %s input.sdf existing.sdf picked.sdf number_to_pick"%sys.argv[0])
        print("Selected diverse number_to_pick compounds from input on top of compounds already picked")
    else:
        input_file = sys.argv[1]
        existing_file = sys.argv[2]
        picked_file = sys.argv[3]
        num_to_pick = int(sys.argv[4])

        fp_list = []
        seed_mol_list = []
        sd_reader = Chem.SDMolSupplier(existing_file)
        idx = 0
        seed_list = []
        mol_list = []
        for mol in sd_reader:
            mol_name = mol.GetProp("Molecule name")
            seed_mol_list.append(mol_name)
            fp_list.append(GetMorganFingerprint(mol,3))
            seed_list.append(idx)
            mol_list.append(mol)
            idx += 1

        sd_reader = Chem.SDMolSupplier(input_file)
        mol: Chem.rdchem.Mol

        for mol in sd_reader:
            if mol.GetProp("Molecule name") not in seed_mol_list:
                mol_list.append(mol)
                fp_list.append(GetMorganFingerprint(mol,3))

        nfp = len(fp_list)
        print(nfp,len(seed_list))
        def dist_ij(i, j, fps=fp_list):
            return 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])

        picker = MaxMinPicker()
        picked_indices = picker.LazyPick(dist_ij,nfp,num_to_pick,firstPicks=seed_list)
        picks = [mol_list[x] for x in picked_indices]

        sd_writer = Chem.SDWriter(picked_file)
        for mol in picks:
            if mol.GetProp("Molecule name") not in seed_mol_list:
                sd_writer.write(mol)
        sd_writer.close()

