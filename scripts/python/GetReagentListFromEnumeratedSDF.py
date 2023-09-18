#!/usr/bin/env python

from openeye.oechem import *
import sys,os
reagentDirectory = "/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/building_blocks/Filtered_50PSA_200MW/ScrubbedList/"
reagentList = [
"WuXi_Acids_Scrubbed_544.sdf",
"WuXi_Aldehydes_Scrubbed_174.sdf",
"WuXi_Amines_Scrubbed_478.sdf",
"WuXi_BoronicAcidsAndEster_Scrubbed_80.sdf",
"WuXi_SulfonylChloride_Scrubbed_5.sdf",
"WuXi_Halides_Scrubbed_184.sdf"
]

if __name__=="__main__":
    if len(sys.argv)!=2:
        print "Usage:%s product.sdf"%sys.argv[0]
    else:
        if not os.path.isfile(sys.argv[1]):
            print "%s is not there."%sys.argv[1]
            sys.exit(1)

        R1_list = []
        R2_list = []
        ifs = oemolistream()
        ifs.open(sys.argv[1])
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            reagents = mol.GetTitle().split("_")
            if reagents[0] not in R1_list:
                R1_list.append(reagents[0])
            if reagents[1] not in R2_list:
                R2_list.append(reagents[1])
        ifs.close()

        reagentDict = {}
        molDict = {}
        for fname in reagentList:
            reagentType = os.path.basename(fname).split(".")[0].split("_")[1]
            fname = os.path.join(reagentDirectory,fname)
            ifs = oemolistream()
            ifs.open(fname)
            mol = OEGraphMol()
            while OEReadMolecule(ifs,mol):
                molDict[mol.GetTitle()] = OEGraphMol(mol)
                if mol.GetTitle() in R1_list:
                    key = "%s_R1"%reagentType
                    if reagentDict.has_key(key):
                        if mol.GetTitle() not in reagentDict[key]:
                            reagentDict[key].append(mol.GetTitle())
                    else:
                        reagentDict[key] = []
                        reagentDict[key].append(mol.GetTitle())
                else:
                    if mol.GetTitle() in R2_list:
                        key = "%s_R2"%reagentType
                        if reagentDict.has_key(key):
                            if mol.GetTitle() not in reagentDict[key]:
                                reagentDict[key].append(mol.GetTitle())
                        else:
                            reagentDict[key] = []
                            reagentDict[key].append(mol.GetTitle())

        for key in reagentDict.keys():
            ofs = oemolostream()
            ofs.open("%s.sdf"%key)
            for molname in reagentDict[key]:
                mol = molDict[molname]
                OEWriteMolecule(ofs,mol)
            ofs.close()





