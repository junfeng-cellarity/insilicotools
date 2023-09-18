#!/usr/bin/env python
import sys
from openeye.oechem import *
import numpy as np
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
from progressbar import ProgressBar
import csv

if __name__ == "__main__":
    if len(sys.argv)!=3:
        print ("Generate ROC curve for csv file, first row contains column name, ")
        print ("First column is actual data, then prediction 1, prediction 2 etc.")
        print ("%s input.csv cutoff"%sys.argv[0])
    else:
        inputCsv = open(sys.argv[1],"r")
        cutoff = float(sys.argv[2])
        reader = csv.DictReader(inputCsv)
        actual_column = reader.fieldnames[0]
        predictions = reader.fieldnames[1:]
        actual_values = []
        prediction_values = []
        for prediction in predictions:
            prediction_values.append([])
        for row in reader:
            actual_values.append(float(row[actual_column])>cutoff)
            for id,prediction in enumerate(predictions):
                prediction_values[id].append(float(row[prediction]))
        print(prediction_values)
        plt.figure()
        for id,prediction in enumerate(predictions):
            fpr,tpr,_ = roc_curve(actual_values,prediction_values[id])
            roc_auc = auc(fpr,tpr)
            plt.plot(fpr, tpr, label='%s ROC curve (area = %0.2f)' %(prediction,roc_auc))
        plt.plot([0, 1], [0, 1], 'k--')
        plt.xlim([-0.1, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        plt.title('Receiver operating characteristic curve')
        plt.legend(loc="lower right")
        plt.show()


