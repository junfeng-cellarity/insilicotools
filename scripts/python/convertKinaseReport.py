#!/usr/bin/env python
import shlex
import psycopg2
from openeye.oechem import *

def convertMolFileToMol(molString):
    ifs = oemolistream()
    ifs.SetFormat(OEFormat_MDL)
    ifs.openstring(molString)
    mol = OEGraphMol()
    r = OEReadMolecule(ifs,mol)
    ifs.close()
    if r:
        return mol
    else:
        return None

class DBServer:
    def __init__(self):
        self.conn = psycopg2.connect(database='biogen',user='medchem',host='javelin',password='medchem')
        self.cursor = self.conn.cursor()

    def getMolFromDb(self,bio_number):
        self.cursor.execute("select molfile from chemicals where bio_number = '%s'" % bio_number)
        row = self.cursor.fetchone()
        if row is not None:
            return convertMolFileToMol(row[0])
        else:
            return None

    def __del__(self):
        self.cursor.close()
        self.conn.close()

db = DBServer()
tags = []
data_dict ={}
data_tags = []
# tagDict = {}
molList = []
f = open("./data/report.csv","r")
count = 0
hitCountDict = {}
molDict = {}
for line in f:
    if count != 0:
        # data = line.strip().split(",")
        myshlex = shlex.shlex(line.strip())
        myshlex.whitespace = ','
        myshlex.whitespace_split = True
        myshlex.commenters = []
        data = []
        for a in myshlex:
            data.append(a)
        bio_number = data[0].rsplit("-",1)[0]
        kinase_name = data[1]
        concentration = data[4]
        tag_name = "%s_%s"%(data[1],data[4])
        tag_value = data[3]
        if kinase_name not in data_tags:
            data_tags.append(kinase_name)
            hitCountDict[kinase_name] = 0
        if float(tag_value) < 30:
            hitCountDict[kinase_name] += 1
        data_key = "%s_%s"%(bio_number,tag_name)
        data_dict[data_key] = tag_value

        mol = db.getMolFromDb(bio_number)
        if mol is None:
            print "Failed to locate ",
            print bio_number
        else:
            if bio_number not in molList:
                molList.append(bio_number)
            molDict[bio_number] = mol
    else:
        for tag in line.split(","):
            if len(tag) > 0:
                tags.append(tag)
        # for tag in tags:
        #     tagDict[tag] = tags.index(tag)
    count += 1
f.close()

sorted(data_tags)
ofs = oemolostream()
ofs.SetFormat(OEFormat_SDF)
ofs.open("kinase3.sdf")


for bio_number in molList:
    if not molDict.has_key(bio_number):
        continue

    mol = molDict[bio_number]
    numHits = 0
    hasProlem = 0
    dose_response = 0
    hasActivity = 0
    for kinase_name in data_tags:

        concentration1 = 100
        data_key_1 = "%s_%s_%d"%(bio_number,kinase_name,concentration1)
        value1 = data_dict[data_key_1]

        concentration2 = 1000
        data_key_2 = "%s_%s_%d"%(bio_number,kinase_name,concentration2)
        value2 = data_dict[data_key_2]

        v1 = float(value1)
        v2 = float(value2)
        if v1 < 30 or v2 < 30:
            if v1 < 30 and v2 > 30:
                hasProlem = 1
            else:
                OESetSDData(mol,"%s_%s"%(kinase_name,concentration1),value1)
                OESetSDData(mol,"%s_%s"%(kinase_name,concentration2),value2)
                if v1<30 and v2<30:
                    dose_response = 1
                numHits += 1
                hasActivity = 1
    if hasActivity:
        OESetSDData(mol,"Number of Hits","%d"%numHits)
        OESetSDData(mol,"Has Problem","%d"%hasProlem)
        OEWriteMolecule(ofs,mol)
ofs.close()
