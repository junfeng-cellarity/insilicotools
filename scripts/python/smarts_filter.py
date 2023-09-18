#!/usr/bin/env python
__author__ = 'jfeng1'
import sys
from openeye.oechem import *

good_smartsList = [
              ("pyrazole","c1cn[nH1]c1"),
              ("AroC~A~R6N","[#6;aH1]~[*;a]:1:[#7;aH0;X2]:[*;a]:[*;a]:[*;a]:[*;a]:1"),
              ("AroC~R5N","[#6;aH1]~[*;a]:1:[#7;aH0;X2]:[*;a]:[*;a]:[*;a]:1"),
              ("N_O~R6N","[#7,#8;H1]~[*;a]:1:[#7;aH0;X2]:[*;a]:[*;a]:[*;a]:[*;a]:1"),
              ("N_O~R5N","[#7,#8;H1]~[*;a]:1:[#7;aH0;X2]:[*;a]:[*;a]:[*;a]:1"),
                ("AroN~CH2NorO","[#7,#8]-[#6H2]-[*;a]:1:[#7;aH0]:[*;a]:[*;a]:[*;a]:[*;a]:1")
              ]
bad_smartsList = [
    ("Acid","[#6](=O)(O-)"),
    ("Acid","[#6](=O)[OH1]")
]

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Usage: %s input.sdf pass.sdf fail.sdf"%sys.argv[0]
    else:
        bad_ss = []
        for name,smart in bad_smartsList:
            s = OESubSearch(smart)
            bad_ss.append((name,s))

        good_ss = []
        for name,smart in good_smartsList:
            s = OESubSearch(smart)
            good_ss.append((name,s))

        passfile = oemolostream()
        passfile.SetFormat(OEFormat_SDF)
        passfile.open(sys.argv[2])

        failfile = oemolostream()
        failfile.SetFormat(OEFormat_SDF)
        failfile.open(sys.argv[3])

        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            OEAddExplicitHydrogens(mol)
            failed = False
            for name,s in bad_ss:
                OEPrepareSearch(mol,s)
                if s.SingleMatch(mol):
                    failed = True
                    OESetSDData(mol, "Fail_Reason", name)
                    OESuppressHydrogens(mol)
                    OEWriteMolecule(failfile,mol)
                    break
            if not failed:
                failed = True
                for name,s in good_ss:
                    if s.SingleMatch(mol):
                        failed = False
                        OESetSDData(mol,"Pass_Reason",name)
                        break
                if failed:
                    OESuppressHydrogens(mol)
                    OEWriteMolecule(failfile,mol)
                else:
                    OESuppressHydrogens(mol)
                    OEWriteMolecule(passfile,mol)
        ifs.close()
        passfile.close()
        failfile.close()
