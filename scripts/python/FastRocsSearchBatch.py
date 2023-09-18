#!/usr/bin/env python

from openeye.oechem import *
from openeye.oemolprop import *
import sys, os, xmlrpclib, glob, random,operator
from openeye.oegraphsim import *
from sqlitedict import SqliteDict
from progressbar import ProgressBar
import random


FASTROCS_SERVER = "http://seabass.biogen.com:%d/RPC2" % 9528

input_dir = "/Users/jfeng1/Datasets/Parkin/"

# > <Abase_ID>
# CDF775328
#
# > <Abase_batch>
# 5
#
# > <dtag_nr>
# DT2012-0356212
#
# > <plate_ID>
# IBS_384_7362
#
# > <well_ref>
# D001
#
# > <ShapeTanimoto>
# 0.8112
#
# > <ColorTanimoto>
# 0.4912
#
# > <TanimotoCombo>
# 1.3024


def convertMolToMolString(mol):
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_SDF)
    ofs.openstring()
    OEWriteMolecule(ofs,mol)
    molString = ofs.GetString()
    ofs.close()
    return molString


class FastRocsServer:
    def __init__(self):
        self.xmlrpc_server = xmlrpclib.ServerProxy(FASTROCS_SERVER)
        pass

    def search(self,inputSdf,numHits):
        input = xmlrpclib.Binary(inputSdf)
        idx = self.xmlrpc_server.SubmitQuery(input,numHits,"sdf","sdf")
        while True:
            blocking = True
            try:
                current,total = self.xmlrpc_server.QueryStatus(idx,blocking)
            except xmlrpclib.Fault,e:
                print >> sys.stderr,str(e)
                return None
            if total == 0:
                continue
            # print "%i/%i"%(current,total)
            if total <= current:
                break

        result = self.xmlrpc_server.QueryResults(idx)
        return result.data

def searchAllCompounds():
    input_files = ["rocs_queries.sdf"]
    output_files = ["output_biib.sdf"]
    smiDict = {}
    server = FastRocsServer()
    for idx,f in enumerate(input_files):
        input = os.path.join(input_dir, f)
        output = os.path.join(input_dir,output_files[idx])
        ofs = oemolostream()
        ofs.open(output)
        ifs = oemolistream()
        ifs.open(input)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            molString = convertMolToMolString(mol)
            result = server.search(molString,100)
            ifs1 = oemolistream()
            ifs1.SetFormat(OEFormat_SDF)
            ifs1.openstring(result)
            mol1 = OEGraphMol()
            while OEReadMolecule(ifs1,mol1):
                smiles = OEMolToSmiles(mol1)
                if smiDict.has_key(smiles):
                    continue
                else:
                    smiDict[smiles] = 1
                score = float(OEGetSDData(mol1,"TanimotoCombo"))
                if score >= 1.2:
                    oechem.OESetSDData(mol1,"Reference",mol.GetTitle())
                    OEWriteMolecule(ofs,mol1)
        ofs.close()



if __name__ == "__main__":
    searchAllCompounds()

