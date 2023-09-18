import re
__author__ = 'jfeng1'

class Compound:
    def __init__(self, compoundId):
        self.compoundId = compoundId
        self.kinaseHits = []
        self.kinaseAct = {}

    def addKinase(self,kinaseId, activity):
        if kinaseId not in self.kinaseHits:
            self.kinaseHits.append(kinaseId)
            self.kinaseAct[kinaseId] = activity

    def __eq__(self, other):
        return self.compoundId == other.compoundId



class Activity:
    def __init__(self,value,type,unit):
        self.value = value
        self.type = type
        self.unit = unit

    # def __eq__(self, other):
    #     if self.type == other.type and self.unit == other.unit:
    #         return True
    #     else:
    #         return False

    def __str__(self):
        return "%s %s %s"%(self.type,self.unit,self.value)

    def isActive(self):
        v = None
        try:
            v = float(self.value)
        except:
            v = None

        if v is None:
            return False

        if self.type == "IC50" or self.type == "Ki" or self.type == "Kd":
            if self.unit == "uM":
                if v <= 0.1:
                    return True
                else:
                    return False
            if self.unit == "nM":
                if v <= 100:
                    return True
                else:
                    return False

class Kinase:
    def __init__(self, domId,name,classification,ecNo):
        self.domId = domId
        self.name = name
        self.classification = classification
        self.ecNo = ecNo

protein_txt = "/Users/jfeng1/Datasets/kinaseSafari/ks_protein.txt"

kinase = []
kinaseNameDict = {}
for line in open(protein_txt,"r"):
    dom = line.split("\t")[0]
    name = line.split("\t")[1]
    lvl4_class = line.split("\t")[5]
    ecNo = line.split("\t")[9]
    k = Kinase(dom, name, lvl4_class, ecNo)
    kinaseNameDict[dom] = k
    kinase.append(dom)

activity_txt = "/Users/jfeng1/Datasets/kinaseSafari/ks_bioactivity.txt"
compounds = []
compoundDict = {}
domList = []
for line in open(activity_txt,"r"):
    activityId,domId,name,assayType,compoundId,activityType,relation,standardValue,standardUnit,activityComment,chemblActivityId,chemblAssayId,pubMedId = line.split("\t")
    if assayType != 'B':
        continue
    if relation == '>' or relation == '>=':
        continue
    if len(domId.strip()) == 0:
        continue
    k = kinaseNameDict[domId]
    act = Activity(standardValue,activityType,standardUnit)
    if not act.isActive():
        continue
    if compoundDict.has_key(compoundId):
        compoundDict[compoundId].addKinase(domId,act)
    else:
        c = Compound(compoundId)
        c.addKinase(domId,act)
        compounds.append(c)
        compoundDict[compoundId] = c
    if domId not in domList:
        domList.append(domId)

for c in compounds:
    print c.compoundId,
    for k in c.kinaseHits:
        print kinaseNameDict[k], c.kinaseAct[k]

from openeye.oechem import *

sdf_file = "/Users/jfeng1/Datasets/kinaseSafari/ks_compound.sdf"

ofs = oemolostream()
ofs.open("kinase_hits_all_detail2.sdf")
ofs.SetFormat(OEFormat_SDF)

known_smiles = []
known_smiles_file = "/Users/jfeng1/JavaProjects/chembldb/python/known_inhibitor_list.txt"
ifs1 = oemolistream()
ifs1.open(known_smiles_file)
m = OEGraphMol()
while OEReadMolecule(ifs1,m):
    known_smiles.append(OECreateCanSmiString(m).split()[0])

ifs = oemolistream()
ifs.open(sdf_file)
mol = OEGraphMol()
while OEReadMolecule(ifs,mol):
    s = OECreateCanSmiString(mol).split()[0]
    if s in known_smiles:
        print "Found %s"%(s)
        continue
    compoundId = OEGetSDData(mol,"COMPOUND_ID")
    if compoundDict.has_key(compoundId):
        c = compoundDict[compoundId]
        OESetSDData(mol,"Number of Hits","%d"%len(c.kinaseHits))
        kinaseNames = []
        for k in c.kinaseHits:
            kinaseTag = kinaseNameDict[k].name
            if kinaseTag not in kinaseNames:
                kinaseNames.append(kinaseTag)
        kinaseNames.sort()
        v = " ".join(kinaseNames)
        OESetSDData(mol,"Hitting Kinases",v)
        OEWriteMolecule(ofs,mol)
ofs.close()

print "Number of targets: %d"%len(kinase)
print "Number of target covered: %d"%len(domList)
print domList
# for p in kinases:
#     if len(dict[p]) > 0:
#         print p,dict[p]

