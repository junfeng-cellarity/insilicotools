#!/usr/bin/env python

import os
import subprocess
from openeye.oechem import *
import random
import tempfile

LOGP_TAG = "Analiza_based_LogP"
LOGP_SDF = "/Users/jfeng1/Databases/ADME/irak_logd_train.sdf"
PERCENT_TRAINING = 0.80
NUM_OF_EXPERIMENTS = 20

class MyTmpFile:
    def __init__(self,mySuffix):
        tuple = tempfile.mkstemp(suffix = mySuffix, text = False)
        if tuple is not None and len(tuple)==2:
            self.file_handle = tuple[0]
            self.full_pathname = tuple[1]
        else:
            self.file_handle = None
            self.full_pathname = None
        return

    def getFullPath(self):
        return self.full_pathname

    def getFileHandle(self):
        return self.file_handle

    def __del__(self):
        if self.full_pathname is not None:
            os.remove(self.full_pathname)


if __name__ == "__main__":
    ifs = oemolistream()
    ifs.open(LOGP_SDF)
    oemol = OEGraphMol()
    molList = []
    actDict = {}
    while OEReadMolecule(ifs,oemol):
        logp = float(OEGetSDData(oemol,LOGP_TAG))
        molList.append(OEGraphMol(oemol))
        actDict[oemol.GetTitle()] = logp
    ifs.close()

    stddev = 0.0
    stddev_user = 0.0
    n = 0
    for i in range(0,NUM_OF_EXPERIMENTS):
        idList = range(0,len(molList))
        random.shuffle(idList)
        k = int(PERCENT_TRAINING*len(molList))
        training_set = idList[0:k]
        validation_set = idList[k:len(molList)]
        #Writing Training Set
        trainFile = MyTmpFile(".sdf")
        ofs = oemolostream()
        ofs.open(trainFile.getFullPath())
        for molId in training_set:
            mol = molList[molId]
            OEWriteMolecule(ofs,mol)
        ofs.close()

        #Run training
        #/Applications/marvinbeans/bin/cxtrain logp -i irak_logp_2 -t Analiza_based_LogP -a irak_logd_train.sdf
        p = subprocess.Popen(['/Applications/marvinbeans/bin/cxtrain','logp','-i','train_%d'%i, '-t','%s'%LOGP_TAG, '-a','%s'%trainFile.getFullPath()],stdout=subprocess.PIPE)
        p.communicate()
        print "Training done"

        #Writing Validation Set
        validationFile = MyTmpFile(".sdf")
        molNames = []
        ofs = oemolostream()
        ofs.open(validationFile.getFullPath())
        for molId in validation_set:
            mol = molList[molId]
            OEWriteMolecule(ofs,mol)
            molNames.append(mol.GetTitle())
        ofs.close()
        import time
        time.sleep(5)
        #Predicting Validation Set
        predictedUser = {}
        p = subprocess.Popen(['/Applications/marvinbeans/bin/cxcalc',validationFile.getFullPath(), 'logp', '--method','user','--trainingid','train_%d'%i], stdout=subprocess.PIPE)
        output = p.communicate()[0]
        lines = output.split('\n')
        for line in lines[1:]:
            if len(line.strip())==0:
                continue
            id,predicted = line.split('\t')
            id = int(id)
            predictedUser[molNames[id-1]] = predicted

        predictedOriginal = {}
        p = subprocess.Popen(['/Applications/marvinbeans/bin/cxcalc',validationFile.getFullPath(), 'logp'], stdout=subprocess.PIPE)
        output = p.communicate()[0]
        lines = output.split('\n')
        for line in lines[1:]:
            if len(line.strip())==0:
                continue
            id,predicted = line.split('\t')
            id = int(id)
            predictedOriginal[molNames[id-1]] = predicted


        for molId in validation_set:
            mol = molList[molId]
            logp = float(OEGetSDData(mol, LOGP_TAG))
            logp_orig = float(predictedOriginal[mol.GetTitle()])
            logp_user = float(predictedUser[mol.GetTitle()])
            print mol.GetTitle(), round(logp,2), round(logp_orig,2) , round(logp_user,2)
            n += 1
            stddev += (logp_orig-logp)**2
            stddev_user += (logp_user-logp)**2
        print
    import math
    print "Standard LogP RMSD:", math.sqrt(stddev/n), "User LogP RMSD:", math.sqrt(stddev_user/n)


