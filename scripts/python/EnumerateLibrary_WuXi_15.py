#!/usr/bin/env python

from openeye.oechem import *
from openeye.oemolprop import *
from openeye.oedepict import *
import sys, os, xmlrpclib, glob, random,operator


reagentDb = {}

smarts_excluded = {
    "[F]": 6,
    "[Br]":1,
    "[N+](=[OX1])[OX1]":0,
    "[CX3]=[CX3]":0,
    "[CX3](=[OX1])[CX3](=[OX1])":0,
    "[CX3](=[OX1])[CX4]([FX1])([FX1])[FX1]":0,
    "[Cl,Br,I]":3,
    "[Br,I]":2,
    "O-O":0,
    "S-S":0,
    "[N]-[N]":0,
    "[N]=[N]":0,
    "[N][N;H2]":0,
    "[CX4;H2,H3;!r][CX4;H2;!r][CX4;H2;!r][CX4;H2,H3;!r]":0,
    "[!#1;!#6;!#7;!#8;!#9;!#35;!#17;!#53;!#16]":0,
    "[N][S;!$(S(=O)(=O)N)]":0,
    "S(~O)(~O)(~O)":0,
    "[N,S][F,Cl,F,I]":0,
    "a-[O;H1]":1,
    "a-[NX3;H2]":1,
    "[A]=[A]-[A]=[A]":0,
    "a-[SX2;H1]":0,
    "A-[Cl,Br,I]":0,
    "C#C-C#C":0,
    "[N;!r]=[C;!r]-Cl":0,
    "*=C(C-O)C=O":0,
    "[N]#[N+]":0,
    "[N,S,O]=C=N":0,
    "S-C#N":0,
    "C=N=N":0,
    "C-N=O":0,
    "N~N~N":0,
    "[C,O]-[CX3;H1]=O":0,
    "[C;!r]-[C;!r]=[N;!r]-[C;!r]":0,
    "C(C#N)C#N":0,
    "[S;!r]=C":0,
    "C(=O)[F,Cl,Br,I]":0,
    "C(=O)C#N":0,
    "O=C-O-C=O":0,
    "a-O-C(=O)A":0,
    "N=C-Cl":0,
    "[OD1&!R]=C1C(=[OD1&!R])C=CC=C1":0,
    "O=[#6]1[#6]:,=[#6][#6](=O)[#6]:,=[#6]1":0,
    "[CX4;r][OX2;H1].[CX4;r][OX2;H1].[CX4;r][OX2;H1].[CX4;r][OX2;H1]":0,
    "[#7,#8,#16;r3]":0,
    "[NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(N-[a])] ":1
}


class MolPropertyCalculator:
    def __init__(self):
        self.xmlrpc_server = xmlrpclib.Server("http://javelin.corp.biogen.com:9528")
        pass

    def calculateMW(self,mol):
        return OECalculateMolecularWeight(mol)

    def calculateXLogP(self,mol):
        return OEGetXLogP(mol)

    def calculateCLogP(self,mol):
        smiles = OEMolToSmiles(mol)
        return self.xmlrpc_server.CLogP(smiles)

    def calculatePSA(self,mol):
        return OEGet2dPSA(mol)

    def calculateAromaticRings(self,mol):
        return OEGetAromaticRingCount(mol)

    def calculateNumberOfRings(self,mol):
        num_components, component_membership = OEDetermineComponents(mol)
        num_rings = mol.NumBonds() - mol.NumAtoms() + num_components
        return num_rings


    def isPassedSmartsChecking(self,mol):
        keys = smarts_excluded.keys()
        subSearch = OESubSearch()
        for key in keys:
            maxMatch = smarts_excluded[key]
            subSearch.Init(key)
            if maxMatch==0:
                if subSearch.SingleMatch(mol):
                    print >>sys.stderr, key, OEMolToSmiles(mol)
                    return False
            else:
                numMatch  = 0
                for ma in subSearch.Match(mol,True):
                    numMatch += 1
                if numMatch > maxMatch:
                    print >>sys.stderr, key, OEMolToSmiles(mol)
                    return False
        return True



class Enumerator:
    def __init__(self, smarts):
        self.libraryGen = OELibraryGen()
        self.libraryGen.Init(smarts)
        self.libraryGen.SetExplicitHydrogens(False)
        self.libraryGen.SetValenceCorrection(True)
        self.libraryGen.SetRemoveUnmappedFragments(True)
        self.numReactants = self.libraryGen.NumReactants()
        print self.numReactants

    def getProduct(self,reactants):
        if len(reactants) != self.numReactants:
            return None
        else:
            for idx,reactant in enumerate(reactants):
                self.libraryGen.SetStartingMaterial(reactant.oemol,idx)
            for product in self.libraryGen.GetProducts():
                return Product(OEGraphMol(product))


class Product:
    def __init__(self, oemol):
        self.oemol = oemol
        OEDeleteEverythingExceptTheFirstLargestComponent(oemol)
        opts = OEPrepareDepictionOptions()
        opts.SetSuppressHydrogens(True)
        opts.SetAddDepictionHydrogens(True)
        opts.SetClearCoords(True)
        opts.SetPerceiveBondStereo(True)
        OEPrepareDepiction(self.oemol,opts)

        self.smiles = OEMolToSmiles(oemol)
        self.name = oemol.GetTitle()
        calculator = MolPropertyCalculator()
        self.nAroRings = calculator.calculateAromaticRings(oemol)
        self.mw = calculator.calculateMW(oemol)
        self.psa = calculator.calculatePSA(oemol)
        self.r1_name = oemol.GetTitle().split("_")[0]
        self.r2_name = oemol.GetTitle().split("_")[1]
        self.xlogp = calculator.calculateXLogP(oemol)
        self.numOfRings = calculator.calculateNumberOfRings(oemol)
        self.passed = calculator.isPassedSmartsChecking(oemol)
        OEClearSDData(self.oemol)
        OESetSDData(self.oemol,"R1_Name",self.r1_name)
        OESetSDData(self.oemol,"R1_Smiles",reagentDb[self.r1_name].smiles)
        OESetSDData(self.oemol,"R2_Name",self.r2_name)
        OESetSDData(self.oemol,"R2_Smiles",reagentDb[self.r2_name].smiles)
        OESetSDData(self.oemol,"MW","%5.2f"%self.mw)
        OESetSDData(self.oemol,"PSA","%5.2f"%self.psa)
        OESetSDData(self.oemol,"XLogP","%5.2f"%self.xlogp)


    def __eq__(self, other):
        return self.smiles == other.smiles

    def __repr__(self):
        return "%s %s",self.smiles,self.name

    def isValid(self):
        if self.mw<=450 and self.psa<=100 and self.xlogp<=4.0 and self.nAroRings>=1 and self.nAroRings<=3 and self.passed and self.numOfRings<=4:
            return True
        return False

class Reagent:
    def __init__(self, oemol):
        OEDeleteEverythingExceptTheFirstLargestComponent(oemol)
        self.oemol = oemol
        self.name = oemol.GetTitle()
        self.smiles = OEMolToSmiles(oemol)
        self.clusterId = int(OEGetSDData(oemol,"KMeanCluster_MACCS166"))
        self.mw = float(OEGetSDData(oemol,"R_MW"))
        self.psa = float(OEGetSDData(oemol,"R_MW"))
        self.numAroRings = OEGetAromaticRingCount(oemol)
        self.hasAromaticRing = self.numAroRings>=1
        reagentDb[self.name] = self

    def __eq__(self, other):
        return self.smiles == other.smiles

    def __repr__(self):
        return "%s %s",self.smiles,self.name

def test():
    calculator = MolPropertyCalculator()
    enumerator = Enumerator("[#6:1]-[#17,#35,#53].[#6:2]-[#7;X3;H2,H1;!$(NC=O):3]>>[#6:2]-[#7:3]-[#6](=O)-[#6]-1=[#6]-2-[#6]-[#7](-[#6:1])-[#6]-[#6]-[#6]-[#7]-2-[#7]=[#7]-1")

    mol1 = OEGraphMol()
    OESmilesToMol(mol1,"CCCCBr")
    mol2 = OEGraphMol()
    OESmilesToMol(mol2,"C1=CC2=CC=CC=C2C=C1N")
    product = enumerator.getProduct([mol1, mol2])
    print OEMolToSmiles(product)

def getDiverseMolecules(reagents,exitingReagents,numReagents):
    idxList = []
    selectedReagents = []
    selectedClusters = {}
    for idx,reagent in enumerate(reagents):
        idxList.append(idx)
    random.shuffle(idxList)
    while len(selectedReagents) < numReagents:
        for idx in idxList:
            reagent = reagents[idx]
            if reagent not in exitingReagents and reagent not in selectedReagents:
                clusterId = reagent.clusterId
                if not selectedClusters.has_key(clusterId):
                    selectedReagents.append(reagent)
                    selectedClusters[clusterId] = 1
                    if len(selectedReagents) == numReagents:
                        break
        selectedClusters.clear()
    return selectedReagents

def readReagents(reagentFile):
    ifs = oemolistream()
    ifs.open(reagentFile)
    reagents = []
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        reagent = Reagent(OEGraphMol(mol))
        reagents.append(reagent)
    ifs.close()
    return reagents

if __name__ == "__main__":
    reagentDirectory = "/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/building_blocks/Filtered_50PSA_200MW/ScrubbedList/"
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX101344/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX100147/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX100037S1/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX102166/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX105109/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX115003/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX140027/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX141022/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX141651/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX105792/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX105148/*.sma")
    # files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/FirstBatch/WX141311/*.sma")

    reagentDict = {
        "Acid":os.path.join(reagentDirectory,"WuXi_Acids_Scrubbed_544_Clustered.sdf"),
        "Aldehyde":os.path.join(reagentDirectory,"WuXi_Aldehydes_Scrubbed_174_Clustered.sdf"),
        "Amine":os.path.join(reagentDirectory,"WuXi_Amines_Scrubbed_478_Clustered.sdf"),
        "BoronicAcidAndEster":os.path.join(reagentDirectory,"WuXi_BoronicAcidsAndEster_Scrubbed_80_Clustered.sdf"),
        "SulfonylChloride":os.path.join(reagentDirectory,"WuXi_SulfonylChloride_Scrubbed_5.sdf"),
        "Halide":os.path.join(reagentDirectory,"WuXi_Halides_Scrubbed_184_Clustered.sdf")
    }
    reagentListDb = {}
    prefixes = []
    for f in files:
        prefix = os.path.basename(f).split(".")[0].split("_")[1]
        if prefix not in prefixes:
            prefixes.append(prefix)
    print "Prefixes:",prefixes
    numFiles = len(prefixes)
    for f in files:
        prefix = os.path.basename(f).split(".")[0]
        coreName,R1_ReagentType, R2_ReagentType = prefix.split("_")
        reagentListKey = "%s_%s"%(coreName,R1_ReagentType)
        print coreName,R1_ReagentType,R2_ReagentType
        R1_reagents = readReagents(reagentDict[R1_ReagentType])
        R2_reagents = readReagents(reagentDict[R2_ReagentType])
        smarts = open(f,"r").read().strip()
        enumerator = Enumerator(smarts)

        if reagentListDb.has_key(reagentListKey):
            productDict = {}
            selectedReagents_R1 = reagentListDb[reagentListKey]
            selectedReagents_R2 = getDiverseMolecules(R2_reagents,[],50)
            successCountR2 = {}
            for r in selectedReagents_R2:
                successCountR2[r.name] = 0

            for r1 in selectedReagents_R1:
                r1=reagentDb[r1[0]]
                for r2 in selectedReagents_R2:
                    r12 = [r1,r2]
                    product = enumerator.getProduct(r12)
                    if product is not None:
                        productDict["%s_%s"%(r1.name,r2.name)] = product
                        if product.isValid():
                            successCountR2[r2.name] += 1

            goodR2 = []
            for key in successCountR2.keys():
                goodR2.append([key,successCountR2[key]])

            goodR2 = sorted(goodR2,key=operator.itemgetter(1), reverse=True)
            nProducts = 0
            directory = os.path.dirname(f)
            output_sdfname = os.path.join(directory,"%s.sdf"%prefix)
            ofs = oemolostream()
            ofs.open(output_sdfname)
            for r1 in selectedReagents_R1:
                for r2 in goodR2[0:30]:
                    if r2[1] < 5:
                        continue
                    key = "%s_%s" % (r1[0], r2[0])
                    if productDict.has_key(key):
                        product = productDict[key]
                        if product.isValid():
                            OEWriteMolecule(ofs,product.oemol)
                            nProducts += 1
            ofs.close()
            print nProducts
        else:
            selectedReagents_R1 = getDiverseMolecules(R1_reagents,[],100)
            successCountR1 = {}
            for r in selectedReagents_R1:
                successCountR1[r.name] = 0
            print "Selected R1:",len(selectedReagents_R1)
            productDict = {}

            selectedReagents_R2 = getDiverseMolecules(R2_reagents,[],50)
            successCountR2 = {}
            for r in selectedReagents_R2:
                successCountR2[r.name] = 0

            print "Selected R2:",len(selectedReagents_R2)
            for r1 in selectedReagents_R1:
                for r2 in selectedReagents_R2:
                    # print "Enumerating %s and %s "%(r1.name,r2.name),
                    r12 = [r1,r2]
                    product = enumerator.getProduct(r12)
                    if product is not None:
                        productDict["%s_%s"%(r1.name,r2.name)] = product
                        if product.isValid():
                            successCountR1[r1.name] += 1
                            successCountR2[r2.name] += 1

            goodR1 = []
            for key in successCountR1.keys():
                goodR1.append([key,successCountR1[key]])

            goodR2 = []
            for key in successCountR2.keys():
                goodR2.append([key,successCountR2[key]])

            goodR1 = sorted(goodR1,key=operator.itemgetter(1), reverse=True)
            goodR2 = sorted(goodR2,key=operator.itemgetter(1), reverse=True)
            nProducts = 0
            directory = os.path.dirname(f)
            output_sdfname = os.path.join(directory,"%s.sdf"%prefix)
            ofs = oemolostream()
            ofs.open(output_sdfname)
            nR1 = 15
            if numFiles == 1:
                nR1 = 25
            reagentListDb[reagentListKey] = goodR1[0:nR1]
            products = []
            for r1 in goodR1[0:nR1]:
                for r2 in goodR2[0:30]:
                    if r2[1] < 5:
                        continue
                    key = "%s_%s" % (r1[0], r2[0])
                    if productDict.has_key(key):
                        product = productDict[key]
                        if product.isValid() and product not in products:
                            OEWriteMolecule(ofs,product.oemol)
                            products.append(product)
                            nProducts += 1
            ofs.close()
            print nProducts
