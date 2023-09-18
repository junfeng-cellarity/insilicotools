#!/usr/bin/env python
import os,sys
from rdkit import Chem

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage:%s input.sdf output.sdf"%(sys.argv[0]))
        print("add core names to sdf")
        sys.exit(1)
    core_name_sdf = "/Users/junfeng/Datasets/DiversitySelection/CDD_full_corenames_current.sdf"
    core_sd_reader = Chem.SDMolSupplier(core_name_sdf)
    dict = {}
    for mol in core_sd_reader:
        if mol.HasProp("Molecule Name"):
            name = mol.GetProp("Molecule Name")
            if mol.HasProp("Core Name Synonym"):
                core_name = mol.GetProp("Core Name Synonym")
                dict[name] = core_name
            else:
                print("None")

    input_sdf = sys.argv[1]
    sd_reader = Chem.SDMolSupplier(input_sdf)
    output_sdf = sys.argv[2]
    sd_writer = Chem.SDWriter(output_sdf)
    for mol in sd_reader:
        name = mol.GetProp("_Name")
        if name in dict:
            mol.SetProp("Core Name",dict[name])
        sd_writer.write(mol)
    sd_writer.close()

