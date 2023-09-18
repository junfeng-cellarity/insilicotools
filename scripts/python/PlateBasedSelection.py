#!/usr/bin/env python

#from openeye.oechem import *
#from openeye.oemolprop import *
import sys, os, glob, random,operator
#from openeye.oegraphsim import *
from sqlitedict import SqliteDict
from progressbar import ProgressBar
import random
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

'''
def convertMolToMolString(mol):
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_SDF)
    ofs.openstring()
    OEWriteMolecule(ofs,mol)
    molString = ofs.GetString()
    ofs.close()
    return molString
'''

def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    print(l,n)
    for i in range(0, len(l), n):
        yield l[i:i+n]

def getAvgSimilarityScore(fpDict, plateId1, plateId2):
    fpList1 = fpDict[plateId1]
    fpList2 = fpDict[plateId2]
    totalSims = 0.0
    n = 0
    for fp1 in fpList1:
        for fp2 in fpList2:
            n += 1
            totalSims +=  DataStructs.TanimotoSimilarity(fp1,fp2)
    return totalSims/n

def genPlateSimilarityMatrix(idx, plateList, allPlates, result_dict, plateDict):
    print(idx," job started.")
    for plate1 in plateList:
        for plate2 in allPlates:
            score = getAvgSimilarityScore(plateDict,plate1,plate2)
            result_dict["%s-%s"%(plate1,plate2)] = score
    print (idx," job finished.")

def getDiversePlates(plateDb, compoundsDict, nCompounds):
    dict = {}
    db = SqliteDict(plateDb)
    for key in db.keys():
        dict[key] = db[key]
    print ("Database loaded into memory.")
    plates = list(compoundsDict.keys())
    print(plates)
    print("OK")
    random.shuffle(plates)

    selectedPlates = []
    nSelected = 0
    while nSelected<nCompounds:
        print ("No.Plates selected:",len(selectedPlates))
        print ("No. Cmpds selected:", nSelected)
        if len(selectedPlates)==0:
            selectedPlates.append(plates[0])
            nSelected = compoundsDict[plates[0]]
            plates.remove(plates[0])
        else:
            minSumOfSim = 9999999
            selectedPlate = None
            for candidate in plates:
                sumOfSim = 0.0
                for plate in selectedPlates:
                    key = "%s-%s"%(candidate,plate)
                    sumOfSim += dict[key]
                if sumOfSim < minSumOfSim:
                    minSumOfSim = sumOfSim
                    selectedPlate = candidate
            selectedPlates.append(selectedPlate)
            print (selectedPlate)
            plates.remove(selectedPlate)
            nSelected += compoundsDict[selectedPlate]

    for plate in selectedPlates:
        print (plate,compoundsDict[plate])

    print (len(selectedPlates)," plates selected.")
    print (nSelected, " compounds selected.")

    return selectedPlates

def generatePlateSimDb(sdfFile,plateIdCol):
    plateDict = {}
    ifs = Chem.SDMolSupplier(sdfFile)

    # ifs = oemolistream()
    # ifs.open(sdfFile)
    # mol = OEGraphMol()
    # while OEReadMolecule(ifs,mol):
    for mol in ifs:
        plate_id = mol.GetProp(plateIdCol)
        fp = Chem.AllChem.GetMorganFingerprint(mol,3)
        # plate_id = OEGetSDData(mol,plateIdCol)
        # fp = OEFingerPrint()
        # OEMakeFP(fp, mol, OEFPType_Circular)
        if not plate_id in plateDict:
            plateDict[plate_id] = []
        plateDict[plate_id].append(fp)
    #ifs.close()
    plates = sorted(plateDict.keys())
    n_processors = 4
    trunk_size = int(len(plates)/n_processors)
    trunks = list(chunks(plates,trunk_size))
    # import pprint
    # pprint.pprint(trunks)
    import multiprocessing
    manager = multiprocessing.Manager()
    result_dict = manager.dict()
    jobs = []
    for idx,trunk in enumerate(trunks):
         p = multiprocessing.Process(target=genPlateSimilarityMatrix, args=(idx,trunk,plates, result_dict, plateDict))
         jobs.append(p)
         p.start()
    for p in jobs:
        p.join()

    db = SqliteDict("PlateSim.db")
    for plateId in result_dict.keys():
        db[plateId] = result_dict[plateId]
    db.commit()


if __name__ == "__main__":

    #########################Generate inter plates similarity matrix ###############################
    if len(sys.argv)==3:
        generatePlateSimDb(sys.argv[1],sys.argv[2])
    #########################Generate inter plates similarity matrix ###############################
    elif len(sys.argv)==5:
        compoundsDict = {}
        molList = []
        ifs = Chem.SDMolSupplier(sys.argv[1])
        for mol in ifs:
            plate_id = mol.GetProp(sys.argv[3])
            if plate_id in compoundsDict:
                compoundsDict[plate_id] = compoundsDict[plate_id]+1
            else:
                compoundsDict[plate_id] = 1
            molList.append(mol)

        plates = getDiversePlates(sys.argv[2],compoundsDict, 5000)

        ofs = Chem.SDWriter(sys.argv[4])
        for mol in molList:
            plate_id = mol.GetProp(sys.argv[3])
            if plate_id in plates:
                ofs.write(mol)
        for plate in plates:
            print(plate)
    else:
        print("Usage:%s input.sdf PlateIdCol" % sys.argv[0])
        print("Usage:%s input.sdf PlateSim.db plateIdCol output.sdf" % sys.argv[0])


