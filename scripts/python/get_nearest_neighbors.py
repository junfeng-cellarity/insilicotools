#!/usr/bin/env python

import rdkit
import os, sys
from rdkit import Chem
from rdkit.SimDivFilters.SimilarityPickers import TopNOverallPicker

db_sdf = "/Users/junfeng/Database/ESC_ChemDiv150k.sdf"
probe_sdf = "/Users/junfeng/Downloads/ESC_hits.sdf"

if __name__ == "__main__":
    probe_names = set()
    sd_reader_2 = Chem.SDMolSupplier(probe_sdf)
    for mol in sd_reader_2:
        mol_name = mol.GetProp("Molecule Name")
        probe_names.add(mol_name)

    sd_reader = Chem.SDMolSupplier(db_sdf)
    mols = []
    for mol in sd_reader:
        mol_name = mol.GetProp("Molecule Name")
        if mol_name not in probe_names:
            mol.SetProp("_Name", mol_name)
            mols.append(mol)

    probe_fps = []
    for mol in mols:
        fp = Chem.RDKFingerprint(mol)
        fp.id = mol.GetProp("_Name")
        probe_fps.append(fp)

    sd_writer = Chem.SDWriter("/Users/junfeng/ESC_hits_NN.sdf")
    sd_reader_2 = Chem.SDMolSupplier(probe_sdf)
    used_nns = set()
    for mol in sd_reader_2:
        probe_fp = Chem.RDKFingerprint(mol)
        picker = TopNOverallPicker(numToPick=5,probeFps=[probe_fp],dataSet=probe_fps)
        for fp,score in picker:
            id = fp.id
            if id not in used_nns:
                used_nns.add(id)
                mol.SetProp("Nearest Neighbor",id)
                mol.SetProp("Tanimoto","%5.2f"%score)
                sd_writer.write(mol)
                break
    sd_writer.close()
