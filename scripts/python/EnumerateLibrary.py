#!/usr/bin/env python

from openeye.oechem import *
from openeye.oemolprop import *
from openeye.oedepict import *
import sys, os, xmlrpclib, glob, random,operator
from sqlitedict import SqliteDict

reagentDb = {}
reagent_usage_db = SqliteDict("/Users/jfeng1/Databases/ReagentUsage.db",autocommit=True)

# Rotatable bond # < 7
# Number of rings >= 3

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
    "[CX4;!$(C~[N,O,S]);!r][CX4;!$(C~[N,O,S]);!r][CX4;!$(C~[N,O,S]);!r][CX4;!$(C~[N,O,S]);!r]":0,
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
    "[NX3;H2,H1,H0;!$(NC=O);!$(NS=O);!$(N[a])] ":1
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
        return float(self.xmlrpc_server.CLogP(smiles))

    def calculatePSA(self,mol):
        return OEGet2dPSA(mol)

    def calculateAromaticRings(self,mol):
        return OEGetAromaticRingCount(mol)

    def calculateNumberOfRings(self,mol):
        num_components, component_membership = OEDetermineComponents(mol)
        num_rings = mol.NumBonds() - mol.NumAtoms() + num_components
        return num_rings

    def calculateNumberOfRotors(self,mol):
        return OECount(mol, OEIsRotor())



    def isPassedSmartsChecking(self,mol):
        keys = smarts_excluded.keys()
        subSearch = OESubSearch()
        for key in keys:
            maxMatch = smarts_excluded[key]
            subSearch.Init(key)
            if maxMatch==0:
                if subSearch.SingleMatch(mol):
                    print >>sys.stderr, key, OEMolToSmiles(mol)
                    return "%s>0"%key
            else:
                numMatch  = 0
                for ma in subSearch.Match(mol,True):
                    numMatch += 1
                if numMatch > maxMatch:
                    print >>sys.stderr, key, OEMolToSmiles(mol)
                    return "%s>%d"%(key,maxMatch)
        return None



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
        args = oemol.GetTitle().split("_")
        self.nComp = len(args)
        self.r1_name = args[0]
        if self.nComp == 2:
            self.r2_name = args[1]
        # self.clogp = calculator.calculateCLogP(oemol)
        self.xlogp = calculator.calculateXLogP(oemol)
        self.passed = calculator.isPassedSmartsChecking(oemol)
        self.numOfRings = calculator.calculateNumberOfRings(oemol)
        OEClearSDData(self.oemol)
        OESetSDData(self.oemol,"R1_Name",self.r1_name)
        OESetSDData(self.oemol,"R1_Smiles",reagentDb[self.r1_name])
        if self.nComp == 2:
            OESetSDData(self.oemol,"R2_Name",self.r2_name)
            OESetSDData(self.oemol,"R2_Smiles",reagentDb[self.r2_name])
        OESetSDData(self.oemol,"MW","%5.2f"%self.mw)
        OESetSDData(self.oemol,"PSA","%5.2f"%self.psa)
        OESetSDData(self.oemol,"XLogP","%5.2f"%self.xlogp)

    def __eq__(self, other):
        return self.smiles == other.smiles

    def __repr__(self):
        return "%s %s",self.smiles,self.name

    def isValid(self):
        valid = False
        if self.mw<=450 and self.psa<=100 and self.xlogp<=4.0 and self.nAroRings>=1 and self.nAroRings<=3 and self.numOfRings>=3 and self.numOfRings<=4 and self.passed is None:
            valid = True
        if not valid:
            reasons = []
            if self.mw>450:
                reasons.append("MW>450")
            if self.psa>100:
                reasons.append("PSA>100")
            if self.xlogp>4:
                reasons.append("XLogP>4")
            if self.nAroRings==0:
                reasons.append("NumAromaticRings=0")
            if self.nAroRings>3:
                reasons.append("NumAromaticRings>3")
            if self.numOfRings>4:
                reasons.append("NumRings>4")
            if self.passed is not None:
                reasons.append(self.passed)
            OESetSDData(self.oemol,"FailedReason"," ".join(reasons))
        return valid

class Reagent:
    def __init__(self, oemol):
        OEDeleteEverythingExceptTheFirstLargestComponent(oemol)
        self.oemol = oemol
        self.name = oemol.GetTitle()
        self.smiles = OEMolToSmiles(oemol)
        self.clusterId = int(OEGetSDData(oemol,"KMeanCluster_MACCS166"))
        self.mw = float(OEGetSDData(oemol,"R_MW"))
        self.psa = float(OEGetSDData(oemol,"R_PSA"))
        self.numAroRings = OEGetAromaticRingCount(oemol)
        self.hasAromaticRing = self.numAroRings>=1
        reagentDb[self.name] = self.smiles

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

def test2():
    mol = OEGraphMol()
    OESmilesToMol(mol,"C(C(CC)C)N1CC2C(C1)C(NC21CCN(CC1)C=1C=CC(C(F)(F)F)=CC1)=O")
    calculator = MolPropertyCalculator()
    print calculator.calculateNumberOfRings(mol)


def getDiverseMolecules(reagents,exitingReagents,numReagents):
    idxList = []
    selectedReagents = []
    selectedClusters = {}
    for idx,reagent in enumerate(reagents):
        idxList.append(idx)
    random.shuffle(idxList)
    while len(selectedReagents) < numReagents:
        print len(selectedReagents)
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

def readReagents(reagentFile,reagentDict):
    ifs = oemolistream()
    ifs.open(reagentFile)
    reagents = []
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        if reagent_usage_db.has_key(mol.GetTitle()) and reagent_usage_db[mol.GetTitle()]>10:
            continue
        reagent = Reagent(OEGraphMol(mol))
        reagents.append(reagent)
        if not reagentDict.has_key(reagent.name):
            reagentDict[reagent.name] = reagent
    ifs.close()
    return reagents



if __name__ == "__main__":
    if True:

        # reagentDirectory = "/Users/jfeng1/WuXi_Library/WuXi_Library_15Cores/building_blocks/Filtered_50PSA_200MW/ScrubbedList/";
        reagentDirectory = "/Users/jfeng1/WuXi_Library/BuildingBlocks";
        # files = glob.glob("/Users/jfeng1/Datasets/Enamine/First10_Enumeration/*.sma")
        # files = glob.glob("/Users/jfeng1/Datasets/Enamine/First10_Enumeration/*.sma")
        #files = glob.glob("/Users/jfeng1/WuXi_Library/WuXi_95Cores/All95Cores/*.sma")
        files = glob.glob("/Users/jfeng1/WuXi_Library/NewEnumeration/*.sma")
        # reagentFileDict = {
        #     "Acid":os.path.join(reagentDirectory,"WuXi_Acids_Scrubbed_544_Clustered.sdf"),
        #     "Aldehyde":os.path.join(reagentDirectory,"WuXi_Aldehydes_Scrubbed_174_Clustered.sdf"),
        #     "Amine":os.path.join(reagentDirectory,"WuXi_Amines_Scrubbed_478_Clustered.sdf"),
        #     "BoronicAcidAndEster":os.path.join(reagentDirectory,"WuXi_BoronicAcidsAndEster_Scrubbed_80_Clustered.sdf"),
        #     "SulfonylChloride":os.path.join(reagentDirectory,"WuXi_SulfonylChloride_Scrubbed_5.sdf"),
        #     "Halide":os.path.join(reagentDirectory,"WuXi_Halides_Scrubbed_184_Clustered.sdf")
        # }
        # Acid_484.sdf
        # Aldehyde_118.sdf
        # Amine_335.sdf
        # BoronicAcidAndEster_75.sdf
        # Halide_186.sdf
        # SulfonylChloride_59.sdf

        reagentFileDict = {
            "Acid":os.path.join(reagentDirectory,"Acid_484.sdf"),
            "Aldehyde":os.path.join(reagentDirectory,"Aldehyde_118.sdf"),
            "Amine":os.path.join(reagentDirectory,"Amine_335.sdf"),
            "BoronicAcidAndEster":os.path.join(reagentDirectory,"BoronicAcidAndEster_75.sdf"),
            "SulfonylChloride":os.path.join(reagentDirectory,"SulfonylChloride_59.sdf"),
            "Halide":os.path.join(reagentDirectory,"Halide_186.sdf")
        }

        maxDict = {
            "Acid":484,
            "Aldehyde":118,
            "Amine":335,
            "BoronicAcidAndEster":75,
            "SulfonylChloride":59,
            "Halide":186
        }
        reagentListDb = {}
        reagentDict = {}
        for reagentType in reagentFileDict.keys():
            reagentListDb[reagentType] = readReagents(reagentFileDict[reagentType], reagentDict)

        r1Dict = {}
        r2Dict = {}
        coreList = set()
        smirksDict = {}
        directory = None
        for f in files:
            if directory is None:
                directory = os.path.dirname(f)
            smirks = open(f,"r").read().strip()
            prefix = os.path.basename(f).split(".")[0]
            args = prefix.split("_")
            if len(args) == 2:
                coreName,R1_ReagentType = args
                smirksDict["%s_%s"%(coreName,R1_ReagentType)] = smirks
            else:
                coreName,R1_ReagentType, R2_ReagentType = args
                smirksDict["%s_%s_%s"%(coreName,R1_ReagentType,R2_ReagentType)] = smirks
                if not r2Dict.has_key(coreName):
                    r2Dict[coreName] = []
                if R2_ReagentType not in r2Dict[coreName]:
                    r2Dict[coreName].append(R2_ReagentType)
                    r2Dict[coreName].sort()
            # coreList.append(coreName)
            coreList.add(coreName)

            if not r1Dict.has_key(coreName):
                r1Dict[coreName] = []
            if R1_ReagentType not in r1Dict[coreName]:
                r1Dict[coreName].append(R1_ReagentType)
                r1Dict[coreName].sort()

        # print r1Dict
        # print r2Dict
        # print coreList

        compound_dict = []
        for coreName in coreList:
            print coreName
            r1List = r1Dict[coreName]
            if r2Dict.has_key(coreName):
                r2List = r2Dict[coreName]
            else:
                r2List = None

            for R1_ReagentType in r1List:
                if r2List is None:
                    productDict = {}
                    R1_reagents = reagentListDb[R1_ReagentType]
                    maxR1 = max(5,int(0.7 * maxDict[R1_ReagentType]))
                    # maxR1 = 30
                    selectedReagents_R1 = getDiverseMolecules(R1_reagents, [], maxR1)
                    # print "Selected R1:",len(selectedReagents_R1)
                    successCountR1 = {}
                    for r in R1_reagents:
                        successCountR1[r.name] = 0
                    for r1 in selectedReagents_R1:
                        enumerator = Enumerator(smirksDict["%s_%s"%(coreName,R1_ReagentType)])
                        # maxR2 = min(30,maxDict[R2_ReagentType])
                        r11 = [r1]
                        product = enumerator.getProduct(r11)
                        if product is not None:
                            productDict[r1.name] = product
                        else:
                            print r1.name," is not working."
                            continue
                        if product.isValid():
                            successCountR1[r1.name] += 1

                    goodR1 = []
                    for key in successCountR1.keys():
                        goodR1.append([key,successCountR1[key]])

                    goodR1 = sorted(goodR1,key=operator.itemgetter(1), reverse=True)
                    print goodR1
                    output_sdfname = os.path.join(directory,"%s_%s.sdf"%(coreName,R1_ReagentType))
                    ofs = oemolostream()
                    ofs.open(output_sdfname)
                    products = []
                    numR1 = 20
                    clusters_r1 = []
                    for r1 in goodR1[0:numR1]:
                        reagent = reagentDict[r1[0]]
                        if reagent.clusterId in clusters_r1:
                            continue
                        else:
                            clusters_r1.append(reagent.clusterId)
                        productName = r1[0]
                        if productDict.has_key(productName):
                            product = productDict[productName]
                            if product is not None and product not in products and product not in compound_dict:
                                if product.isValid():
                                    OESetSDData(product.oemol,"Valid","YES")
                                else:
                                    OESetSDData(product.oemol,"Valid","No")
                                OEWriteMolecule(ofs,product.oemol)
                                products.append(product)
                                compound_dict.append(product)
                            else:
                                print "failed...."
                    ofs.close()
                    continue
                else:
                    productDict = {}
                    R1_reagents = reagentListDb[R1_ReagentType]
                    maxR1 = max(5,int(0.7 * maxDict[R1_ReagentType]))
                    # maxR1 = 30
                    selectedReagents_R1 = getDiverseMolecules(R1_reagents, [], maxR1)
                    # print "Selected R1:",len(selectedReagents_R1)
                    successCountR1 = {}
                    successCountR2 = {}
                    for r in R1_reagents:
                        successCountR1[r.name] = 0
                    for R2_ReagentType in r2List:
                        R2_reagents = reagentListDb[R2_ReagentType]
                        for r in R2_reagents:
                            successCountR2[r.name] = 0

                    for r1 in selectedReagents_R1:
                        for R2_ReagentType in r2List:
                            enumerator = Enumerator(smirksDict["%s_%s_%s"%(coreName,R1_ReagentType,R2_ReagentType)])
                            maxR2 = max(5,int(0.7 * maxDict[R2_ReagentType]))
                            # maxR2 = min(30,maxDict[R2_ReagentType])
                            R2_reagents = reagentListDb[R2_ReagentType]
                            selectedReagents_R2 = getDiverseMolecules(R2_reagents,[],maxR2)
                            for r2 in selectedReagents_R2:
                                r12 = [r1,r2]
                                product = enumerator.getProduct(r12)
                                if product is not None:
                                    productDict["%s_%s"%(r1.name,r2.name)] = product
                                else:
                                    print r1.name,r2.name," is not working."
                                    continue
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
                    print goodR1
                    print goodR2
                    output_sdfname = os.path.join(directory,"%s_%s_%s.sdf"%(coreName,R1_ReagentType,"_".join(r2List)))
                    ofs = oemolostream()
                    ofs.open(output_sdfname)
                    products = []
                    numR1 = 8
                    numR2 = 10
                    clusters_r1 = []
                    for r1 in goodR1:
                        print r1[0], successCountR1[r1[0]]
                        reagent = reagentDict[r1[0]]
                        if reagent.clusterId in clusters_r1:
                            continue
                        else:
                            clusters_r1.append(reagent.clusterId)

                        clusters_r2 = []
                        for r2 in goodR2:
                            reagent = reagentDict[r2[0]]
                            if reagent.clusterId in clusters_r2:
                                continue
                            else:
                                clusters_r2.append(reagent.clusterId)
                            productName = "%s_%s" % (r1[0], r2[0])
                            if productDict.has_key(productName):
                                product = productDict[productName]
                                if product is not None and product not in products and product not in compound_dict:
                                    if product.isValid():
                                        OESetSDData(product.oemol,"Valid","YES")
                                    else:
                                        OESetSDData(product.oemol,"Valid","No")
                                    OEWriteMolecule(ofs,product.oemol)
                                    products.append(product)
                                    compound_dict.append(product)
                                else:
                                    print "failed...."
                    ofs.close()






