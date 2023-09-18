#!/usr/bin/env python
from openeye.oechem import *
from openeye.oemolprop import *
import sys, os, xmlrpclib, glob, random,operator
import json

class MolPropertyCalculator:
    def __init__(self):
        self.xmlrpc_server = xmlrpclib.Server("http://javelin.corp.biogen.com:9528")
        pass

    def calculateMW(self,mol):
        mw = OECalculateMolecularWeight(mol)
        OESetSDData(mol,"MW", "%5.3f"%mw)

    def calculateCLogP(self,molList):
        smiList = []
        for molId,mol in enumerate(molList):
            smiles = OEMolToSmiles(mol)
            smiList.append("%s %d"%(smiles,molId))
        input = "\n".join(smiList)
        result = self.xmlrpc_server.CLogPBatch(input)
        dict = json.loads(result)
        for molId,mol in enumerate(molList):
            molKey = "%d" % molId
            if dict.has_key(molKey):
                OESetSDData(mol,"CLogP","%s"%dict[molKey])

    def calculatePSA(self,mol):
        psa = OEGet2dPSA(mol)
        OESetSDData(mol,"PSA","%5.3f"%psa)

    def calculateCFp3(self,mol):
        fp3 = OEGetFractionCsp3(mol)
        OESetSDData(mol,"fp3","%5.3f"%fp3)

    def calculateAromaticRings(self,mol):
        nRing = OEGetAromaticRingCount(mol)
        OESetSDData(mol,"nAromaticRings","%3d"%nRing)

    def calculateNumberOfRings(self,mol):
        num_components, component_membership = OEDetermineComponents(mol)
        num_rings = mol.NumBonds() - mol.NumAtoms() + num_components
        OESetSDData(mol,"nRings","%3d"%num_rings)

    def calculateNumberOfRotors(self,mol):
        nRotors = OECount(mol, OEIsRotor())
        OESetSDData(mol,"nRotors","%3d"%nRotors)

    def calculateHBD_HBA(self,mol):
        OEThrow.SetLevel(OEErrorLevel_Warning)
        hba_name = "hydrogen-bond acceptors"
        hbd_name = "hydrogen-bond donors"
        oefilter = OEFilter(OEFilterType_Lead)
        ostr = oeosstream()
        pwnd = False
        oefilter.SetTable(ostr, pwnd)
        headers = ostr.str().split(b'\t')
        ostr.clear() # remove this row from the stream
        oefilter(mol)
        fields = ostr.str().decode("UTF-8").split('\t')
        tmpdct = dict(zip(headers, fields))
        hba = tmpdct[hba_name]
        hbd = tmpdct[hbd_name]
        OESetSDData(mol,hba_name,hba)
        OESetSDData(mol,hbd_name,hbd)

if __name__=="__main__":
    if len(sys.argv)<3:
        print "Usage:%s input.sdf output.sdf"%sys.argv[0]
        print "Generating common descriptors for molecules"
    else:
        calculator = MolPropertyCalculator()
        inputSdf = sys.argv[1]
        outputSdf = sys.argv[2]
        molList = []
        ifs = oemolistream()
        ifs.open(inputSdf)
        mol = OEGraphMol()
        ofs = oemolostream()
        ofs.open(sys.argv[2])
        while OEReadMolecule(ifs,mol):
            newMol = OEGraphMol(mol)
            calculator.calculatePSA(newMol)
            calculator.calculateNumberOfRings(newMol)
            calculator.calculateMW(newMol)
            calculator.calculateAromaticRings(newMol)
            calculator.calculateCFp3(newMol)
            calculator.calculateNumberOfRotors(newMol)
            calculator.calculateHBD_HBA(newMol)
            molList.append(newMol)
        calculator.calculateCLogP(molList)
        for mol in molList:
            OEWriteMolecule(ofs,mol)
        ifs.close()
        ofs.close()