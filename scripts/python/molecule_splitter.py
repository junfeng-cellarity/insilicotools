#!/usr/bin/env python
import os, sys
from openeye.oechem import *

def get_ring_nbr(nbr1: OEAtomBase, nbr2:OEAtomBase):
    nbr: OEAtomBase
    if nbr2.IsInRing():
        return 1
    count = 0
    for nbr in nbr2.GetAtoms():
        if nbr == nbr1:
            continue
        if nbr.IsInRing():
            count += 1
        else:
            count += get_ring_nbr(nbr2,nbr)
    return count

def is_atom_linker(atom1: OEAtomBase):
    if atom1.IsInRing():
        return False
    else:
        count = 0
        for nbr_atm in atom1.GetAtoms():
            count += get_ring_nbr(atom1,nbr_atm)
        if count == 2:
            return True
        else:
            return False

if __name__ == "__main__":
    ifs = oemolistream()
    ifs.open("/Users/junfeng/Pharmaron_HiBit_110K.sdf")
    #ifs.open("/Users/junfeng/error_splitted.sdf")
    row_dicts = {}
    mol = OEGraphMol()
    mol_list = []
    while OEReadMolecule(ifs,mol):
        orginal_mol = OEGraphMol(mol)
        atom: OEAtomBase
        linker_anchor_atoms = []
        for atom in mol.GetAtoms():
            if is_atom_linker(atom):
                for nbr in atom.GetAtoms():
                    if nbr.IsInRing():
                        linker_anchor_atoms.append(atom)

        linker_bonds = []
        bond: OEBondBase
        for bond in mol.GetBonds():
            if bond.IsInRing():
                continue
            else:
                if bond.GetBgn().IsInRing() and bond.GetEnd().IsInRing():
                    linker_bonds.append(bond)

        broken_bond_atoms = []
        for bond in linker_bonds:
            new_atm1 = mol.NewAtom(OEElemNo_Au)
            new_atm2 = mol.NewAtom(OEElemNo_Au)
            new_bond1 = mol.NewBond(new_atm1, bond.GetBgn())
            new_bond2 = mol.NewBond(new_atm2, bond.GetEnd())
            mol.DeleteBond(bond)
            broken_bond_atoms.append((bond.GetBgnIdx(),bond.GetEndIdx()))

        for atom in linker_anchor_atoms:
            for bond in atom.GetBonds():
                if bond.GetBgn() == atom and bond.GetEnd().IsInRing():
                    new_atm1 = mol.NewAtom(OEElemNo_Au)
                    new_atm2 = mol.NewAtom(OEElemNo_Au)
                    new_bond1 = mol.NewBond(new_atm1, bond.GetBgn())
                    new_bond2 = mol.NewBond(new_atm2, bond.GetEnd())
                    mol.DeleteBond(bond)
                    broken_bond_atoms.append((bond.GetBgnIdx(), bond.GetEndIdx()))
                elif bond.GetEnd() == atom and bond.GetBgn().IsInRing():
                    new_atm1 = mol.NewAtom(OEElemNo_Au)
                    new_atm2 = mol.NewAtom(OEElemNo_Au)
                    new_bond1 = mol.NewBond(new_atm1, bond.GetBgn())
                    new_bond2 = mol.NewBond(new_atm2, bond.GetEnd())
                    mol.DeleteBond(bond)
                    broken_bond_atoms.append((bond.GetBgnIdx(), bond.GetEndIdx()))
        OEAssignMDLHydrogens(mol)
        OEGenerate2DCoordinates(mol)
        OEAssignAromaticFlags(mol)

        numparts, partlist = OEDetermineComponents(mol)
        pred = OEPartPredAtom(partlist)

        parts = []
        for i in range(1, numparts + 1):
            pred.SelectPart(i)
            partmol = OEGraphMol()
            OESubsetMol(partmol, mol, pred)
            parts.append(partmol)

        basic_amine_smarts = "[NX3;!$(NC=O);!$(NS=O);!$(Na);!$(NC=[*]);!$(N[!#6;!#1;!#79])]"
        basic_amine_pred = OEMatchAtom(basic_amine_smarts)
        au_pred = OEMatchAtom("[Au]")
        ring_pred = OEMatchAtom("[r]")
        order_dict = {}

        for id,partmol in enumerate(parts):
            basic_N_count = OECount(partmol,basic_amine_pred)
            Au_count = OECount(partmol,au_pred)
            ring_count = OECount(partmol,ring_pred)

            # print(OEMolToSmiles(partmol),basic_N_count,Au_count,ring_count)
            if Au_count == 2 and ring_count > 0:
                if 0 not in order_dict:
                    order_dict[0] = []
                order_dict[0].append(id)
            if Au_count == 1 and ring_count > 0 and basic_N_count == 0:
                order_dict[1] = id
            if Au_count >= 1 and basic_N_count > 0:
                order_dict[2] = id

        keys = sorted(order_dict.keys())
        sorted_parts = []
        for key in keys:
            if key==0:
                sub_mol = OEGraphMol()
                for idx in order_dict[key]:
                    OEAddMols(sub_mol,parts[idx])
                sorted_parts.append(sub_mol)
            else:
                sorted_parts.append(parts[order_dict[key]])
        row_dicts[mol.GetTitle()] = sorted_parts
        mol_list.append(orginal_mol)

    import csv
    keys = [0,1,2]
    field_names = ["SMILES","Name"]
    for id in keys:
        field_names.append("R%d_SMILES"%(id+1))
    field_names_act = ["Name","Act"]
    with open("/Users/junfeng/rgroup.csv","w") as csv_file, open("/Users/junfeng/rgroup_act.csv","w") as act_file:
        act_writer = csv.DictWriter(act_file,fieldnames=field_names_act)
        act_writer.writeheader()
        csv_writer = csv.DictWriter(csv_file,fieldnames=field_names)
        csv_writer.writeheader()
        for id,mol in enumerate(mol_list):
            act_row_dict = {}
            act_row_dict["Name"] = mol.GetTitle()
            act_row_dict["Act"] = float(OEGetSDData(mol,"Pharmaron_HiBiT_Normalozed"))
            #act_row_dict["Act"] = float(OEGetSDData(mol,"Act"))
            act_writer.writerow(act_row_dict)

            row_dict = {}
            row_dict["SMILES"] = OEMolToSmiles(mol)
            row_dict["Name"] = mol.GetTitle()
            for id in keys:
                parts = row_dicts[mol.GetTitle()]
                if id< len(parts):
                    row_dict["R%d_SMILES"%(id+1)] = OEMolToSmiles(row_dicts[mol.GetTitle()][id])
            csv_writer.writerow(row_dict)
    csv_file.close()
    act_file.close()

    # ofs = oemolostream()
    # ofs.open("/Users/junfeng/splitted.sdf")
    # for partmol in sorted_parts:
    #     OEWriteMolecule(ofs,partmol)
    # ofs.close()