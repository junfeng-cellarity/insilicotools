__author__ = 'jfeng1'
from openeye.oechem import *

f = open("/Users/jfeng1/224.list","r")
myList = []
for line in f:
    myList.append(line.strip())

ifs = oemolistream()
ifs.open("/cluster/database/databases/Evotec/Registered/evotec_biib_set.sdf")
mol = OEGraphMol()

while OEReadMolecule(ifs,mol):
    bio_number = OEGetSDData(mol,"BIO Number")
    if bio_number in myList:
        lot_number = OEGetSDData(mol,"Lot number")
        print "%s-%s"%(bio_number,lot_number)
        myList.remove(bio_number)


