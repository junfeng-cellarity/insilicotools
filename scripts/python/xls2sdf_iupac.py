#!/usr/bin/env python
from openpyxl import load_workbook
import glob
import os,sys
from openeye.oechem import *
import psycopg2
import psycopg2.extras
import traceback
from parse import parse
from openeye.oeiupac import *

if __name__ == "__main__":
    if len(sys.argv)!=4:
        print "Usage:%s input.xls output.sdf columnId(first column 0)"%sys.argv[0]
        print "First row contains data tag names"
        print "BIO-NUMBER is assumed to be the first column"
        sys.exit(0)

    OEThrow.SetLevel(OEErrorLevel_Error)
    xlsFile = sys.argv[1]
    sdfFile = sys.argv[2]
    columnId = int(sys.argv[3])

    ofs = oemolostream()
    ofs.open(sdfFile)
    dict = {}
    wb = load_workbook(xlsFile)
    active_sheet = wb.get_active_sheet()
    colNames = []
    first_row = active_sheet.rows[0]
    for c in first_row:
        colName = c.value
        if colName is not None:
            colNames.append(str(colName))
    for row in active_sheet.rows[1:]:
        mol = OEGraphMol()
        for colId,colName in enumerate(colNames):
            OESetSDData(mol,colName,str(row[colId].value))
            if colId == columnId:
                iupac_name = row[colId].value
                if iupac_name is None or len(iupac_name) == 0:
                    continue
                OEParseIUPACName(mol, iupac_name)
                mol.SetTitle(iupac_name)
        OEWriteMolecule(ofs,mol)

    ofs.close()

