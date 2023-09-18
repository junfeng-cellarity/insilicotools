#!/usr/bin/env python
from openpyxl import load_workbook
import glob
import os,sys
from openeye.oechem import *
import psycopg2
import psycopg2.extras
import traceback
from parse import parse

class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='biogen2',user='medchem',host='javelin',password='medchem')

    def getMolFromBio(self, bio_number):
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select m from rdk.mols where bio_number = %s",(bio_number,))
            result = cursor.fetchone()
            if result is not None:
                smiles = result['m']
                mol = OEGraphMol()
                ifs = oemolistream()
                ifs.SetFormat(OEFormat_SMI)
                ifs.openstring(smiles)
                OEReadMolecule(ifs,mol)
                return mol

            else:
                return None
        except:
            cursor.execute("rollback")
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()



if __name__ == "__main__":
    if len(sys.argv)!=3:
        print "Usage:%s input.xls output.sdf"%sys.argv[0]
        print "First row contains data tag names"
        print "Smiles must be one of the column names"
        sys.exit(0)

    xlsFile = sys.argv[1]
    sdfFile = sys.argv[2]

    ofs = oemolostream()
    ofs.open(sdfFile)

    wb = load_workbook(xlsFile)
    active_sheet = wb.get_active_sheet()
    colNames = []
    first_row = active_sheet.rows[0]
    smiColId = -1
    for id,c in enumerate(first_row):
	print id,c.value
        colName = c.value
        if colName == None:
            colName = "Col%d"%id
        if colName.upper() == "SMILES":
            smiColId = id
        if colName is not None:
            colNames.append(str(colName))

    for row in active_sheet.rows[1:]:
        mol = OEGraphMol()
        OESmilesToMol(mol,str(row[smiColId].value))
        for colId,colName in enumerate(colNames):
            # OESetSDData(mol,colNames[colId],row[colId].value.decode("utf-8"))
            OESetSDData(mol,colNames[colId],str(row[colId].value))
        OEWriteMolecule(ofs,mol)

    ofs.close()

