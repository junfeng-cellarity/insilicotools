#!/usr/bin/env python
import sys,os
from openeye.oechem import *
from openpyxl import load_workbook
KEY = "Alias_ID"
KEY2 = "Alias Id"
columnsInSdf = ["Tube_barcode","Tube_Net_Amount_mg","Tube_location","Plate_ID"]
columnsInXls = ["container ID","quantity","container position",	"Rack ID"]

if __name__ == "__main__":
    if len(sys.argv)!=4:
        print "Usage:%s fromDonna.xlsx vendor.sdf output.xlsx"%sys.argv[0]
    else:
        xlsFile = sys.argv[1]
        sdfFile = sys.argv[2]
        outputXlsFile = sys.argv[3]

        molDict = {}
        ifs = oemolistream()
        ifs.open(sdfFile)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            aliasId = OEGetSDData(mol,KEY)
            molDict[aliasId] = OEGraphMol(mol)

        colNameDict = {}
        wb = load_workbook(xlsFile)
        sheet = wb.get_active_sheet()
        for rowid,row in enumerate(sheet.rows):
            if rowid ==0:
                for colId,col in enumerate(row):
                    colNameDict[colId] = col.value
            else:
                rowDict = {}
                for colId,col in enumerate(row):
                    colName = colNameDict[colId]
                    rowDict[colName] = col
                aliasId = rowDict[KEY2].value
                mol = molDict[aliasId]
                for id,col in enumerate(columnsInSdf):
                    rowDict[columnsInXls[id]].value = OEGetSDData(mol,col)
        wb.save(sys.argv[3])






