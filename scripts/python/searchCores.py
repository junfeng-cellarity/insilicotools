#!/usr/bin/env python
from openeye.oechem import *
from openeye.oemedchem import *
import sys
CORE_SDF = "/Users/jfeng1/Datasets/CoresCollection/cores.sdf"

if len(sys.argv) != 3:
    print "Usage: %s input.sdf output.sdf"%sys.argv[0]
    sys.exit(1)

ifs = oemolistream()
ifs.open(CORE_SDF)
coreMol = OEGraphMol()
coreNames = []
coreSmilesDict = {}
coreMolDict = {}
coreDict = {}
result = {}

while OEReadMolecule(ifs,coreMol):
    coreName = coreMol.GetTitle()
    coreNames.append(coreName)
    result[coreName] = []
    subsearch = OESubSearch(coreMol,OEExprOpts_DefaultAtoms,OEExprOpts_DefaultBonds)
    coreDict[coreName] = subsearch
    coreSmilesDict[coreName] = OEMolToSmiles(coreMol)
    coreMolDict[coreName] = OEGraphMol(coreMol)
ifs.close()

ofs = oemolostream()
ofs.open(sys.argv[2])

ifs = oemolistream()
ifs.open(sys.argv[1])
mol = OEGraphMol()
while OEReadMolecule(ifs,mol):
    newmol = OEGraphMol(mol)
    found = False
    for coreName in coreNames:
        print coreName
        subsearch = coreDict[coreName]
        coreMol = coreMolDict[coreName]
        OEPrepareSearch(mol,subsearch)
        if subsearch.SingleMatch(mol):
            OESetSDData(newmol,"MatchedCore",coreName)
            OESetSDData(newmol,"MatchedCoreSmiles",coreSmilesDict[coreName])
            OEWriteMolecule(ofs,newmol)
            result[coreName].append(newmol.GetTitle())
            found = True
            break
        else:
            print "No hits."

        OETheFunctionFormerlyKnownAsStripSalts(coreMol)
        OEFindRingAtomsAndBonds(coreMol)
        for frag in OEGetRingChainFragments(coreMol):
            fragatompred = OEIsAtomMember(frag.GetAtoms())
            fragbondpred = OEIsBondMember(frag.GetBonds())

            fragment = OEGraphMol()
            adjustHCount = True
            OESubsetMol(fragment, coreMol, fragatompred, fragbondpred, adjustHCount)
            smiles = OEMolToSmiles(fragment)
            if smiles.rfind("1") != -1 and smiles.rfind("2") != -1:
                subsearch2 = OESubSearch()
                subsearch2.Init(smiles)
                OEPrepareSearch(mol,subsearch2)
                if subsearch2.SingleMatch(mol):
                    OESetSDData(newmol,"coreSmiles_%s"%coreName,smiles)
                    OESetSDData(newmol,"coreSmilesHits_%s"%coreName,coreSmilesDict[coreName])
                    found = True

        for frag in OEGetBemisMurcko(coreMol):
            fragatompred = OEIsAtomMember(frag.GetAtoms())
            fragbondpred = OEIsBondMember(frag.GetBonds())
            fragment = OEGraphMol()
            adjustHCount = True
            OESubsetMol(fragment, coreMol, fragatompred, fragbondpred, adjustHCount)
            smiles = OEMolToSmiles(fragment)
            for oerole in frag.GetRoles():
                if oerole.GetName() == "Framework":
                    subsearch2 = OESubSearch()
                    subsearch2.Init(smiles)
                    OEPrepareSearch(mol,subsearch2)
                    if subsearch2.SingleMatch(mol):
                        OESetSDData(newmol,"frameworkSmiles_%s"%coreName,smiles)
                        OESetSDData(newmol,"frameworkHitsSmiles_%s"%coreName,coreSmilesDict[coreName])
                        found = True

    if found:
        OEWriteMolecule(ofs,newmol)


ifs.close()
ofs.close()

