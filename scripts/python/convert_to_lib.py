#!/usr/bin/env python

from openeye.oechem import *
import re,sys,os


def get_pseudo_name(name):
    name = name.strip()
    if name.startswith("C"):
        return re.sub(r'[C]', r'Q', name)
    if name.startswith("N"):
        return re.sub(r'[N]', r'QN', name)
    return "Q"


class LibAtom:
    def __init__(self, atom, molecule):
        self.atom = atom
        self.molecule = molecule
        self.x = self.molecule.GetCoords(self.atom)[0]
        self.y = self.molecule.GetCoords(self.atom)[1]
        self.z = self.molecule.GetCoords(self.atom)[2]
        self.nbr_names = []
        # self.nbr_names = []
        if not atom.IsHydrogen():
            for nbr in self.atom.GetAtoms():
                self.nbr_names.append("%-5s" % nbr.GetName().strip())
        else:
            for nbr in self.atom.GetAtoms():
                h_count = 0
                for nbr1 in nbr.GetAtoms():
                    if nbr1.IsHydrogen():
                        h_count = h_count + 1
                if h_count >= 2:
                    q_name = get_pseudo_name(nbr.GetName())
                    self.nbr_names.append("%-5s" % nbr.GetName().strip())
                    self.nbr_names.append("%-5s" % "-")
                    self.nbr_names.append("%-5s" % "-")
                    self.nbr_names.append("%-5s" % "-")
                    self.nbr_names.append("%-5s" % q_name)

    def get_name(self):
        return self.atom.GetName()

    def get_symbol(self):
        return OEGetAtomicSymbol(self.atom.GetAtomicNum())

    def print_atom(self, idx, output):
        #       print("1234 1234512345123451234567890123456789012345678901234567890  12345123451234512345")
        print("%4d %-5s%-5s%5d%10.4f%10.4f%10.4f%10.4f  %s"
              % (idx, self.get_name().strip(), self.get_symbol(),
                 0, 0, self.x, self.y, self.z, "".join(self.nbr_names)),file=output)


def print_pseudoatom(idx, atom_name, output):
    #       9 QN1  PSEUD    0    0.0000    0.0000    0.0000    0.0000
    print("%4d %-5s%-5s%5d%10.4f%10.4f%10.4f%10.4f" % (idx, atom_name, "PSEUD", 0, 0.0, 0.0, 0.0, 0.0),file=output)


def print_atom(lib_atom, atom_dict, atom_list, output):
    if lib_atom in atom_list:
        return
    lib_atom.print_atom(len(atom_list) + 1,output)
    atom_list.append(lib_atom)
    hydrogen_atoms = []
    for nbr in lib_atom.atom.GetAtoms():
        if nbr.IsHydrogen():
            nbr_libatom = atom_dict[nbr.GetIdx()]
            hydrogen_atoms.append(nbr_libatom)
            print_atom(nbr_libatom, atom_dict, atom_list, output)
    if len(hydrogen_atoms) >= 2:
        q_name = get_pseudo_name(lib_atom.get_name()).strip()
        atom_list.append(q_name)
        print_pseudoatom(len(atom_list) + 1, q_name, output)
    for nbr in lib_atom.atom.GetAtoms():
        nbr_libatom = atom_dict[nbr.GetIdx()]
        if not nbr.IsHydrogen():
            print_atom(nbr_libatom, atom_dict, atom_list, output)


if __name__ == "__main__":
    if len(sys.argv)!=3:
        print("Usage:%s input.pdb output.lib"%sys.argv[0])
    else:
        pdb_input = sys.argv[1]
        lib_output = sys.argv[2]
        output = open(lib_output,"w")
        atom_list = []
        ifs = oemolistream()
        ifs.open(pdb_input)
        pdb_mol = OEGraphMol()
        OEReadMolecule(ifs, pdb_mol)
        ifs.close()
        print("RESIDUE   1",file=output)
        atom_dict = {}
        for atom in pdb_mol.GetAtoms():
            libatom = LibAtom(atom, pdb_mol)
            atom_dict[atom.GetIdx()] = libatom

        for atom in pdb_mol.GetAtoms():
            print_atom(atom_dict[atom.GetIdx()], atom_dict, atom_list, output)
