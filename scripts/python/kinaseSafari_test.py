import re
__author__ = 'jfeng1'

kinases = ['AAK1',
           'BARK',
           'CRIK',
           'DMPK',
           'GRK',
           'LATS',
           'MAST',
           'MRCK',
           'NDR',
           'PKG',
           'PRKY',
           'RIOK',
           'RSKL',
           'SBK1',
           'TYRO',
           'ULK',
           'WNK']


class Compound:
    def __init__(self, compoundId):
        self.compoundId = compoundId
        self.kinaseHits = []
        self.kinaseAct = {}

    def addKinase(self,kinaseId, activity):
        if not kinaseId in self.kinaseHits:
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

        if self.type == "IC50" or self.type == "Potency" or self.type == "Ki" or self.type == "Kd":
            if self.unit == "uM":
                if v < 5:
                    return True
                else:
                    return False
            if self.unit == "nM":
                if v < 5000:
                    return True
                else:
                    return False


protein_txt = "/Users/jfeng1/Datasets/kinaseSafari/ks_protein.txt"

kinaseNameDict = {}
dict = {}
for p in kinases:
    dict[p] = []

for line in open(protein_txt,"r"):
    for p in kinases:
        m = re.search(r"(%s)"%p,line)
        if m is not None:
            print m.group()
            dom = line.split("\t")[0]
            name = line.split("\t")[1]
            dict[p].append(dom)
            kinaseNameDict[dom] = name


activity_txt = "/Users/jfeng1/Datasets/kinaseSafari/ks_bioactivity.txt"
compounds = []
compoundDict = {}
for line in open(activity_txt,"r"):
    activityId,domId,name,assayType,compoundId,activityType,relation,standardValue,standardUnit,activityComment,chemblActivityId,chemblAssayId,pubMedId = line.split("\t")
    act= Activity(standardValue,activityType,standardUnit)
    if not act.isActive():
        continue
    for p in kinases:
        for domId1 in dict[p]:
            if domId == domId1:
                if compoundDict.has_key(compoundId):
                    compoundDict[compoundId].addKinase(domId,act)
                else:
                    c = Compound(compoundId)
                    c.addKinase(domId,act)
                    compounds.append(c)
                    compoundDict[compoundId] = c

for c in compounds:
    print c.compoundId,
    for k in c.kinaseHits:
        print kinaseNameDict[k], c.kinaseAct[k]


from openeye.oechem import *

sdf_file = "/Users/jfeng1/Datasets/kinaseSafari/ks_compound.sdf"

ofs = oemolostream()
ofs.open("kinase_hits.sdf")
ofs.SetFormat(OEFormat_SDF)

ifs = oemolistream()
ifs.open(sdf_file)
mol = OEGraphMol()
while OEReadMolecule(ifs,mol):
    compoundId = OEGetSDData(mol,"COMPOUND_ID")
    if compoundDict.has_key(compoundId):
        c = compoundDict[compoundId]
        OESetSDData(mol,"Number of Hits","%d"%len(c.kinaseHits))
        for k in c.kinaseHits:
            kinaseTag = kinaseNameDict[k]
            tag1 = "%s type"%kinaseTag
            tag2 = "%s value"%kinaseTag
            OESetSDData(mol,tag1,"%s %s"%(c.kinaseAct[k].type,c.kinaseAct[k].unit))
            OESetSDData(mol,tag2,"%s"%c.kinaseAct[k].value)
        OEWriteMolecule(ofs,mol)
ofs.close()


# for p in kinases:
#     if len(dict[p]) > 0:
#         print p,dict[p]

