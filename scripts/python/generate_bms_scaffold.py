#!/usr/bin/env python
import rdkit
from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol
import sys
if __name__=="__main__":
    if len(sys.argv) != 3:
        print("Usage: %s input.sdf output.sdf"%sys.argv[0])
        print("Generate Bemis-Murcko Scaffold")
    else:
        input_fname = sys.argv[1]
        output_fname = sys.argv[2]
        sd_reader = Chem.SDMolSupplier(input_fname)
        sd_writer = Chem.SDWriter(output_fname)
        mol: rdkit.Chem.rdchem.Mol
        for mol in sd_reader:
            scaffold_mol = GetScaffoldForMol(mol)
            scaffold_smi = Chem.MolToSmiles(scaffold_mol)
            mol.SetProp("MurckoScaffold",scaffold_smi)
            sd_writer.write(mol)
        sd_writer.close()
