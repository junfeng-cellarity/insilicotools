#!/usr/bin/env python

from openeye.oechem import *

three_letter_dict = {
    "A": ["Ala","Alanine","Hydrophobic"],
    "R": ["Arg","Arginine","Positive Charged"],
    "N": ["Asn","Asparagine","Polar"],
    "D": ["Asp","Aspartic acid","Negative Charged"],
    "C": ["Cys","Cysteine","Polar"],
    "E": ["Glu","Glutamic acid","Negative Charged"],
    "Q": ["Gln","Glutamine","Polar"],
    "G": ["Gly","Glycine","Hydrophobic"],
    "H": ["His","Histidine","Polar"],
    "I": ["Ile","Isoleucine","Hydrophobic"],
    "L": ["Leu","Leucine","Hydrophobic"],
    "K": ["Lys","Lysine","Positive Charged"],
    "M": ["Met","Methionine","Polar"],
    "F": ["Phe","Phenylalanine","Hydrophobic"],
    "P": ["Pro","Proline","Hydrophobic"],
    "S": ["Ser","Serine","Polar"],
    "T": ["Thr","Threonine","Polar"],
    "W": ["Trp","Tryptophan","Polar"],
    "Y": ["Tyr","Tyrosine","Polar"],
    "V": ["Val","Valine","Hydrophobic"]
}

class AminoAcid:
    def __init__(self, one_letter):
        self.one_letter = one_letter
        if three_letter_dict.has_key(one_letter):
            details = three_letter_dict[one_letter]
            self.three_letter = details[0]
            self.name = details[1]
            self.type = details[2]
        else:
            self.three_letter = one_letter
            self.name = one_letter
            self.type = "unknown"

    def __str__(self):
        return self.name


class Peptide:
    def __init__(self, sequence):
        self.sequence = list(sequence)
        self.aa_list = []
        for aa in self.sequence:
            self.aa_list.append(AminoAcid(aa))
        self.name = sequence

    def set_name(self,name):
        self.name = name

    def get_seq(self):
        return "".join(self.sequence)

    def mutate(self,idx,aa):
        other_peptide = Peptide(self.get_seq())
        other_peptide.sequence[idx] = aa
        other_peptide.aa_list[idx] = AminoAcid(aa)
        return other_peptide

    def getAA(self,idx):
        return self.aa_list[idx].one_letter

    def __str__(self):
        return ">Id %s\n"%self.name+"".join(self.sequence)



def generate_mutation():
    ofs = oemolostream()
    ofs.open("/Users/jfeng1/fasta.sdf")
    peptide = Peptide("DAEFRHDSG")
    amino_acids = three_letter_dict.keys()
    id = 0
    print peptide.getAA(3),peptide.getAA(4)
    for aa in amino_acids:
        if aa == peptide.getAA(3):
            continue
        for ab in amino_acids:
            if ab == peptide.getAA(5):
                continue
            new_peptide = peptide.mutate(3,aa).mutate(5,ab)
            mol = OEGraphMol()
            ifs = oemolistream()
            ifs.SetFormat(OEFormat_FASTA)
            ifs.openstring(str(new_peptide))
            OEReadMolecule(ifs,mol)
            OEWriteMolecule(ofs,mol)
    ofs.close()
    #
    # ofs = oemolostream()
    # ofs.openstring()
    # ofs.SetFormat(OEFormat_SDF)
    # OEWriteMolecule(ofs,mol)
    # print ofs.GetString()

if __name__ == "__main__":
    peptide = Peptide("DAEFRHDSG")
    ifs = oemolistream()
    ifs.SetFormat(OEFormat_FASTA)
    ifs.openstring(str(peptide))
    mol = OEGraphMol()
    OEReadMolecule(ifs,mol)
    ofs = oemolostream()
    ofs.open("/Users/jfeng1/peptide_inverted.sdf")
    n = 1
    for atm in mol.GetAtoms():
        # atm = OEAtomBase()
        if atm.IsChiral() and atm.IsCarbon():
            print n,atm.GetName()
            newMol = OEGraphMol(mol)
            new_atm = newMol.GetAtom(OEHasAtomIdx(atm.GetIdx()))
            v = []
            for nbr in new_atm.GetAtoms():
                v.append(nbr)
            stereo = new_atm.GetStereo(v, OEAtomStereo_Tetrahedral)
            if stereo == OEAtomStereo_LeftHanded:
                new_atm.SetStereo(v,OEAtomStereo_Tetra,OEAtomStereo_Right)
            else:
                if stereo == OEAtomStereo_Right:
                    new_atm.SetStereo(v,OEAtomStereo_Tetra,OEAtomStereo_Left)
            newMol.SetTitle("%s_%d"%(mol.GetTitle(),n))
            n+=1
            OEWriteMolecule(ofs,newMol)
    ofs.close()
    ifs.close()

