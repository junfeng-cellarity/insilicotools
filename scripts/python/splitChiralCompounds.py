#!/usr/bin/env python
from openeye.oechem import *
import os
import sys

if __name__ == "__main__":
    if len(sys.argv)!=2:
        print "Usage:%s input.sdf "%(sys.argv[0])
    else:
        ifs = oemolistream()
        ifs.open(sys.argv[1])

        mol = OEGraphMol()
        dict = []
        ofs_nonchiral = oemolostream()
        ofs_nonchiral.open("nonchiral1.sdf")
        ofs_chiral_pure = oemolostream()
        ofs_chiral_pure.open("chiral_pure.sdf")
        ofs_chiral_mix = oemolostream()
        ofs_chiral_mix.open("chiral__mix.sdf")
        ofs_error = oemolostream()
        ofs_error.open("chiral_error.sdf")
        while OEReadMolecule(ifs,mol):
            chiral_flag_on = False
            if OEMDLGetParity(mol):
                chiral_flag_on = True
            molIsChiral = False
            molChiralDefined = True
            chiral_pure = True
            for atom in mol.GetAtoms():
                chiral = atom.IsChiral() and atom.IsCarbon()
                if not chiral:
                    continue
                else:
                    molIsChiral = True
                stereo = OEAtomStereo_Undefined
                if atom.HasStereoSpecified(OEAtomStereo_Tetrahedral):
                    v = []
                    for nbr in atom.GetAtoms():
                        v.append(nbr)
                    stereo = atom.GetStereo(v, OEAtomStereo_Tetrahedral)
                if stereo == OEAtomStereo_Undefined:
                    molChiralDefined = False
                else:
                    if not chiral_flag_on:
                        hasGroup = False
                        for group in mol.GetGroups(OEHasGroupType(OEGroupType_MDLOrStereo)):
                            hasGroup = True
                            if atom not in group.GetAtoms():
                                chiral_pure = False
                                break
                        if not hasGroup:
                            chiral_pure = False
                            break

            if not molIsChiral:
                OEWriteMolecule(ofs_nonchiral,mol)
            else:
                if molChiralDefined:
                    if chiral_flag_on:
                        OEWriteMolecule(ofs_chiral_pure,mol)
                    else:
                        if chiral_pure:
                            OEWriteMolecule(ofs_chiral_pure,mol)
                        else:
                            OEWriteMolecule(ofs_chiral_mix,mol)
                else:
                    if chiral_flag_on:
                        OEWriteMolecule(ofs_error,mol)
                    else:
                        OEWriteMolecule(ofs_chiral_mix,mol)

        ofs_chiral_pure.close()
        ofs_chiral_mix.close()
        ofs_nonchiral.close()
        ofs_error.close()

        ifs.close()
