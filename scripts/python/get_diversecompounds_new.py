#!/usr/bin/env python
import rdkit
from rdkit import Chem,DataStructs
import os,sys
from rdkit.SimDivFilters.rdSimDivPickers import MaxMinPicker
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprint

if __name__ == "__main__":
    if len(sys.argv)!=4:
        print("Usage: %s input.sdf output.sdf picked.sdf"%sys.argv[0])
        print("Selected diverse compounds from input")
    else:
        min_cmpds_per_cluster = 1
        num_cmpds = 5000
        scaffold_set = set()
        scaffold_dict = {}
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        picked_file = sys.argv[3]
        scaffold_tag = "MurckoScaffoldSmiles"
        sd_reader = Chem.SDMolSupplier(input_file)
        mol: Chem.rdchem.Mol
        for mol in sd_reader:
            scaffold_smi = mol.GetProp(scaffold_tag)
            scaffold_set.add(scaffold_smi)
            if scaffold_smi not in scaffold_dict:
                scaffold_dict[scaffold_smi] = set()
            scaffold_dict[scaffold_smi].add(mol)

        mol_list = []
        for scaffold_smiles in scaffold_set:
            mini_mol_list = list(scaffold_dict[scaffold_smiles])
            if len(mini_mol_list)<=min_cmpds_per_cluster:
                for mol in mini_mol_list:
                    mol_list.append(mol)
            else:
                mini_fp_list = []
                def mini_dist_ij(i, j, fps=mini_fp_list):
                    return 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])
                for mol in mini_mol_list:
                    fp = GetMorganFingerprint(mol,3)
                    mini_fp_list.append(fp)
                else:
                    mini_nfp = len(mini_fp_list)
                    mini_picker = MaxMinPicker()
                    mini_picked_indices = mini_picker.LazyPick(mini_dist_ij,mini_nfp,min_cmpds_per_cluster,seed=3423)
                    mini_picked  = [mini_mol_list[x1] for x1 in mini_picked_indices ]
                    for mol in mini_picked:
                        mol_list.append(mol)

        sd_writer = Chem.SDWriter(output_file)
        for mol in mol_list:
            sd_writer.write(mol)
        sd_writer.close()

        fp_list = []
        for mol in mol_list:
            fp_list.append(GetMorganFingerprint(mol,3))
        nfp = len(fp_list)

        def dist_ij(i, j, fps=fp_list):
            return 1 - DataStructs.TanimotoSimilarity(fps[i], fps[j])

        picker = MaxMinPicker()
        picked_indices = picker.LazyPick(dist_ij,nfp,num_cmpds,seed=7455)
        picks = [mol_list[x] for x in picked_indices]

        sd_writer = Chem.SDWriter(picked_file)
        for mol in picks:
            sd_writer.write(mol)
        sd_writer.close()

