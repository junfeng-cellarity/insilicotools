#!/usr/bin/env python
from openeye.oechem import *
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from progressbar import ProgressBar

sdfFile = "/Users/jfeng1/TTBK1/ttbk1_result_glidescore_all.sdf"

db = OEMolDatabase(sdfFile)
nCompounds = db.GetMaxMolIdx()

ifs = oemolistream()
ifs.open(sdfFile)
mol = OEGraphMol()
isActive= []
glideScore = []
n = 0
progress = ProgressBar(maxval=nCompounds).start()
while OEReadMolecule(ifs,mol):
    if OEGetSDData(mol,"has_glide_score")=="YES":
        active = 0
        activity = float(OEGetSDData(mol, "%Activity"))
        if activity <50:
            active = 1
        score = -float(OEGetSDData(mol,"glide_score"))
        isActive.append(active)
        glideScore.append(score)
    progress.update(n)
    n += 1

print (isActive,glideScore)
fpr,tpr,_ = roc_curve(isActive,glideScore)
roc_auc = auc(fpr,tpr)


plt.figure()
plt.plot(fpr, tpr, label='ROC curve (area = %0.2f)' % roc_auc)
plt.plot([0, 1], [0, 1], 'k--')
plt.xlim([0.0, 1.0])
plt.ylim([0.0, 1.05])
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic curve')
plt.legend(loc="lower right")
plt.show()


