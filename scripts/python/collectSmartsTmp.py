#!/usr/bin/env python
from openeye.oechem import *
import traceback

smarts_file = "/Users/jfeng1/faf_structure_alerts.txt"
AlertNames = ["Frequents Hitters","Aggregators","Toxicophores","PAINS"]
StructureAlertDict = {}
read_seq = 0
currentAlert = None
currentName = None
f = open(smarts_file, "r")
for s in f:
    pattern = s.strip()
    if pattern in AlertNames:
        StructureAlertDict[pattern] = []
        currentAlert = pattern
        read_seq = 0
    else:
        if read_seq == 0:
            read_seq += 1
            currentName = pattern
        else:
            try:
                subsearch = OESubSearch()
                result = subsearch.Init(pattern)
                read_seq = 0
                if result:
                    print "%s\t%s\t%s\t"%(pattern,currentName,currentAlert)
            except:
                traceback.print_exc()
                read_seq = 0
                pass

f.close()