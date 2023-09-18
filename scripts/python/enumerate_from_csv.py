#!/usr/bin/env python

#from __future__ import print_function
from openeye.oechem import *
from parse import *
import sys

#Example of CSV file to use as input
#name	R0 smiles	R1 smiles	R2 smiles	R3 smiles	R4 smiles
#VC-0001102130	[H]N(c1cc([*:2])ncn1)[*:1]	[H]c1nn(C([H])([H])C([H])([H])[H])c([H])c1[*:1]	[H]c1c([H])c2c(c([H])c1([*:2]))C([H])([H])N(C1([H])C([H])([H])OC1([H])[H])C([H])([H])C([H])([H])C2([H])[*:3]	[H]N(C(=O)[*:4])[*:3]	[H]C([H])([H])C(OC1([H])C([H])([H])N([*:4])C1([H])[H])(C([H])([H])[H])C([H])([H])[H]


def connect_molecule(smiles_array):
    mol = OEGraphMol()
    for smiles in smiles_array:
        frag = OEGraphMol()
        OEParseSmiles(frag, smiles)
        OEAddMols(mol,frag)
    return mol

if __name__=="__main__":
    if len(sys.argv)!=3:
        print "Usage:%s input.csv output.sdf" % sys.argv[0]
        sys.exit(0)

    titles = []
    reagentDict = {}
    propertyDict = {}
    csv = open(sys.argv[1],"r").read()
    for lineno,line in enumerate(csv.splitlines()):
        if lineno==0:
            titles = line.split(",")
        else:
            fragments = line.split(",")
            for idx,fragment in enumerate(fragments):
                rs = parse("R{} smiles", titles[idx])
                if rs is not None:
                    rg_id = rs[0]
                    if rg_id not in reagentDict:
                        reagentDict[rg_id] = []
                    reagentDict[rg_id].append(fragment)
                else:
                    if titles[idx] not in propertyDict:
                        propertyDict[titles[idx]] = []
                    propertyDict[titles[idx]].append(fragment)
    num_cmpds = len(reagentDict['0'])
    rgroup_ids = reagentDict.keys()
    rgroup_ids.sort()
    ofs = oemolostream()
    ofs.open(sys.argv[2])
    for i in range(0,num_cmpds):
        smi_array = []
        for j in rgroup_ids:
            smi_array.append(reagentDict[j][i])
        mol = connect_molecule(smi_array)
        atom_map = {}
        atom_to_delete = []
        for atom in mol.GetAtoms():
            if atom.GetMapIdx() != 0:
                map_idx = atom.GetMapIdx()
                if map_idx not in atom_map:
                    atom_map[map_idx] = []
                for nbr in atom.GetAtoms():
                    atom_map[map_idx].append(nbr)
                    break
                atom_to_delete.append(atom)
        for atomList in atom_map.values():
            if len(atomList)==2:
                mol.NewBond(atomList[0],atomList[1])
            else:
                print >> sys.stderr, "error in mol %d"%i
        for atom in atom_to_delete:
            mol.DeleteAtom(atom)
        for property in propertyDict.keys():
            OESetSDData(mol,property,propertyDict[property][i])
        OEAssignMDLHydrogens(mol)
        OEWriteMolecule(ofs,mol)
    ofs.close()

