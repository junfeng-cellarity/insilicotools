#!/usr/bin/env python

from openeye.oechem import *
import sys,glob,os
from openpyxl import load_workbook

UPPER_DIRECTORY = "/Users/jfeng1/WuXi_Library/building_blocks/Filtered_50PSA_200MW/"
DIRECTORY = "/Users/jfeng1/WuXi_Library/building_blocks/Filtered_50PSA_200MW/ScrubbedList"
filenames = [
    "WuXi_Acids_652.sdf",
    "WuXi_Aldehyde_236.sdf",
    "WuXi_Amines_531.sdf",
    "WuXi_BoronicAcidsAndEster_81.sdf",
    "WuXi_Halide_205.sdf",
    "WuXi_SulphonyChloride_5.sdf"
]

molDB = {}

for f in filenames:
    ifs = oemolistream()
    ifs.open(os.path.join(UPPER_DIRECTORY,f))
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        molDB[mol.GetTitle()] = OEGraphMol(mol)
    ifs.close()

dir = os.path.join(DIRECTORY,"*.xlsx")

files = glob.glob(dir)

for f in files:
    wb = load_workbook(f)
    progress = 0
    sdf_basename = os.path.basename(f).split(".")[0]
    for sheet_name in wb.get_sheet_names():
        active_sheet = wb.get_sheet_by_name(sheet_name)
        mol_numbers = len(active_sheet.rows)-1
        sdf_basename = "%s_Scrubbed_%d.sdf"%(sdf_basename,mol_numbers)
        print sdf_basename
        sdfname = os.path.join(DIRECTORY,sdf_basename)
        ofs = oemolostream()
        ofs.open(sdfname)
        print sdfname
        for rowid,row in enumerate(active_sheet.rows):
            if rowid >0:
                OEWriteMolecule(ofs,molDB[row[1].value])
        ofs.close()

