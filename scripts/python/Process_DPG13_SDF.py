#!/usr/bin/env python
import os,sys
from openeye.oechem import *
smartsDict = {
    "tetracyclic": ["[#6]-1-[#7]-c2[#6,#7;a]cnc3n[#6,#7;a]c(-[#6;a][#6;a]-1)c23"],
    "macrocyclic": ["*-1-*-*-[*;a]:[*;a]:[*;a]-c2[#6,#7;a]nc3nc[#6,#7;a]c(-*-*-1)c23","*-1-*-*-*-c2[#6,#7;a]cnc3n[#6,#7;a]c(-[*;a]:[*;a]:[*;a]-*-*-1)c23","*-1-*-*-*-*-c2[#6,#7;a]cnc3n[#6,#7;a]c(-[*;a]:[*;a]:[*;a]-*-*-1)c23"],
    #"morpholin": ["[#6;a]-[#7]-c1n[#6,#7;a]c2[#6,#7;a][#6,#7;a][#6,#7;a][#6,#7;a]c2[#6,#7;a]1"],
    "5,6-amino_pyrazole":["*-[#7H1]-[*;a]:1:[*;a]:[*;a]:[*;a]:2:[*;a]:[*;a]:[*;a](-*):[*;a]:2:[*;a]:1"],
    "6,6-amino_pyrazole and others":["[#7;h1](-[#6;a])-c1cnn(-[#6])c1"]
}

wh_smartsDict = {
    "Acrylamide":["[#7]-C(=O)-[#6;h1]=[#6;h2]"],
    "Terminal Substituted Acrylamide":["[#7]-C(=O)-[#6;h1]=C-[*;!h]"],
    "NonTerminal Substituted Acrylamide":["[#7]-C(=O)-[#6;h0]=[#6;h2]"],
    "Halo Acrylamide":["[#7]-C(=O)-C([F,Cl,Br])=C"],
    "Double Substituted Acrylamide":["[#7]-C(=O)-[#6]([#6])=C-[*;!h]"],
    "Alkynyl Amide":["NC(=O)-C#C"],
    "HaloAmide":["NC(=O)-[CX4]-[Cl,F,Br]"],
    "Cyano":["C#N"]
}

smarts_list = ["tetracyclic","macrocyclic","5,6-amino_pyrazole","other amino_pyrazole"]
wh_smarts_list = ["Acrylamide","Terminal Substituted Acrylamide","NonTerminal Substituted Acrylamide",
                  "Halo Acrylamide","Double Substituted Acrylamide","Alkynyl Amide","HaloAmide","Cyano"]
PCE_List = [
    "CY-2002081",
    "CY-2002200",
    "CY-2003145",
    "CY-2010777",
    "CY-2011249",
    "CY-2011803",
    "CY-2012229",
    "CY-2012261",
    "CY-2013472",
    "CY-2002262",
    "CY-2005807",
    "CY-2002188",
    "CY-2011243",
    "CY-2004279",
    "CY-2005746",
    "CY-2011241",
    "CY-2014893",
    "CY-2014894",
    "CY-2015027",
    "CY-2015028"
]

sdf_tags = [
    'ID',
    "PCE",
    "Scaffold",
    "Warhead",
    'JAK3 enzyme (nM)',
    'JAK3 cell (nM)',
    'Shake Flask LogD',
    'BTK 1mM ATP GeoMean IC50 (nM)',
    'Hu Blood Stability t1/2 (min)',
    'Rat Blood Stability t1/2 (min)',
    'Hu Hepatocyte t1/2 (min)',
    'Hu Liver Microsomal Stability t1/2 (min)',
    'Rat Hepatocyte t1/2 (min)',
    'Rat Liver Microsomal Stability t1/2 (min)',
    'Dog Hepatocyte t1/2 (min)',
    'Dog Blood Stability t1/2 (min)',
    'Hu PPB Fu (%)',
    'Rat PPB Fu (%)',
    'MDR1 NIH Papp A-B (10^-6 cm/sec)',
    'MDR1 NIH Papp B-A (10^-6 cm/sec)',
    'MDR1 NIH EFFLUX RATIO',
    'Permeability Assay MDCK P-gp KO Cell Line Papp A to B GEOM_MEAN',
    'Pharmaron hERG Patch Clamp (uM)',
    'Human HEPATOCYTES CL int (ul/min/10e6 cells)',
    'Human HEPATOCYTES CL int scaled (ml/min/Kg)',
    'Rat HEPATOCYTES CL int (ul/min/10e6 cells)',
    'Rat HEPATOCYTES CL int scaled (ml/min/Kg)',
    'Blood Stability Human Remaining % at 120 min',
    'Blood Stability Rat Remaining % at 120 min',
    'Pharmaron GSH t1/2 (min)',
    'Pharmaron GSH % Remaining 240 min',
    'Kinetic Solubility pH 7.4 PBS',
    'Solubility Thermodynamic 7.4 PBS GEOM_MEAN',
    'LATEST_ADMET_DATE',
    'LATEST_BIO_DATE',
    'Stereochemistry',
    'Program Code',
    'MDR1 NIH Average Papp (10^-6 cm/sec)',
    'JAK3 1mM ATP Imax (%)',
    'JAK3 Biochemical Preincub IC50 (nM)',
    'JAK3 Biochemical Preincub Imax (%)',
    'JAK3 IC50 Shift',
    'JAK3 Cellular Imax_obs (%)',
    'BTK 1mM ATP GeoMean Imax (%)',
    'MDR1 NIH RECOVERY A-B (pct)',
    'MDR1 NIH BA RECOVERY B-A (pct)',
    'Hu PPB Recovery (%)',
    'Rat PPB Recovery (%)',
]

if __name__ == "__main__":
    if len(sys.argv) == 3:
        input_sdf = sys.argv[1]
        output_sdf = sys.argv[2]
        ifs = oemolistream()
        ifs.open(input_sdf)
        ofs = oemolostream()
        ofs.open(output_sdf)
        mol = OEGraphMol()
        tags = []
        subsearchDict = {}
        for ssName in smartsDict.keys():
            queryList = smartsDict[ssName]
            subsearchDict[ssName] = []
            for id,query in enumerate(queryList):
                ss = OESubSearch(query)
                subsearchDict[ssName].append(ss)

        wh_subsearchDict = {}
        for ssName in wh_smartsDict.keys():
            queryList = wh_smartsDict[ssName]
            wh_subsearchDict[ssName] = []
            for id,query in enumerate(queryList):
                ss = OESubSearch(query)
                wh_subsearchDict[ssName].append(ss)

        while OEReadMolecule(ifs,mol):
            cy_number = OEGetSDData(mol, "ID")
            if cy_number in PCE_List:
                OESetSDData(mol,"PCE","YES")
            else:
                OESetSDData(mol,"PCE","NO")
            foundMatch = False
            for ssName in smarts_list:
                for ss in subsearchDict[ssName]:
                    OEPrepareSearch(mol, ss)
                    if ss.SingleMatch(mol):
                        OESetSDData(mol,"Scaffold", ssName)
                        foundMatch = True
                        break
                if foundMatch:
                    break
            if not foundMatch:
                OESetSDData(mol,"Scaffold","Other")

            foundMatch = False
            for ssName in wh_smarts_list:
                for ss in wh_subsearchDict[ssName]:
                    OEPrepareSearch(mol,ss)
                    if ss.SingleMatch(mol):
                        OESetSDData(mol,"Warhead",ssName)
                        foundMatch = True
                        break
                if foundMatch:
                    break
            if not foundMatch:
                OESetSDData(mol,"Warhead", "Other/No Warhead")

            dataDict = {}
            non_essential_dataDict = {}

            for dp in OEGetSDDataPairs(mol):
                if dp.GetTag() not in sdf_tags:
                    non_essential_dataDict[dp.GetTag()] = dp.GetValue()

            for tag in sdf_tags:
                if OEHasSDData(mol,tag):
                    dataDict[tag] = OEGetSDData(mol,tag)
                else:
                    dataDict[tag] = ""

            OEClearSDData(mol)

            for tag in sdf_tags:
                OESetSDData(mol,tag,dataDict[tag])

            for key in non_essential_dataDict.keys():
                OESetSDData(mol, key, non_essential_dataDict[key])

            OEWriteMolecule(ofs,mol)