__author__ = 'jfeng1'
from openeye.oechem import *
import re
import sys
if len(sys.argv)!=3:
    print "Usage:%s input.sdf output.sdf"%(sys.argv[0])
    sys.exit(1)

ifs = oemolistream()
ifs.open(sys.argv[1])
solvates = ["Hydrate","Isopropanol","Ethanol","Dimethylformamide","Acetonitrile"]
ofs = oemolostream()
ofs.open(sys.argv[2])
ofs.SetFormat(OEFormat_SDF)
mol = OEGraphMol()
while OEReadMolecule(ifs,mol):
    saltString = OEGetSDData(mol, "Salt")
    OEDeleteSDData(mol,"Salt")
    OEDeleteSDData(mol,"Coefficient")
    if saltString is not None:
        salts = re.split(r'(?<=[a-zA-Z ]),(?<![0-9a-zA-Z])',saltString)
        # OESetSDData(mol,"OriginalSalt",saltString)
        # OESetSDData(mol, "Coefficient","%d"%len(salts))
        saltDict = {}
        solvateDict = {}
        for salt in salts:
            salt = salt.strip()
            if salt in solvates:
                if solvateDict.has_key(salt):
                    solvateDict[salt] = solvateDict[salt] + 1
                else:
                    solvateDict[salt] = 1
            else:
                if saltDict.has_key(salt):
                    saltDict[salt] = saltDict[salt] + 1
                else:
                    saltDict[salt] = 1
        id = 1
        for salt in saltDict.keys():
            if id > 1:
                OESetSDData(mol,"Salt%d"%id, salt)
                OESetSDData(mol,"Salt%d_coeff"%id,"%d"%saltDict[salt])
            else:
                OESetSDData(mol,"Salt", salt)
                OESetSDData(mol,"Salt_coeff","%d"%saltDict[salt])
            id += 1
        id = 1
        for solvate in solvateDict.keys():
            if id > 1:
                OESetSDData(mol,"Solvate%d"%id, solvate)
                OESetSDData(mol,"Solvate%d_coeff"%id,"%d"%solvateDict[solvate])
            else:
                OESetSDData(mol,"Solvate", solvate)
                OESetSDData(mol,"Solvate_coeff","%d"%solvateDict[solvate])
            id += 1
    OEWriteMolecule(ofs,mol)
ofs.close()
ifs.close()

