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
        print "BIO-NUMBER is assumed to be the first column"
        sys.exit(0)

    OEThrow.SetLevel(OEErrorLevel_Error)
    xlsFile = sys.argv[1]
    sdfFile = sys.argv[2]

    ofs = oemolostream()
    ofs.open(sdfFile)
    dict = {}
    wb = load_workbook(xlsFile)
    active_sheet = wb.get_active_sheet()
    colNames = []
    first_row = list(active_sheet.rows)[0]
    for c in first_row:
        colName = c.value
        if colName is not None:
            colNames.append(str(colName))
    db = BiogenDb()
    for row in list(active_sheet.rows)[1:]:
        mol = None
        for colId,colName in enumerate(colNames):
            if colId == 0:
                s = str(row[0].value)
                p = parse("BIO-{bionumber}-{lot}",s)
                if p is None:
                    p = parse("BIO-{bionumber}",s)
                if p is None:
                    print s," is empty"
                else:
                    bio_number = "BIO-%s"%p['bionumber']
                    mol = db.getMolFromBio(bio_number)
                if mol is None:
                    print bio_number, " failed to find."
                    mol = OEGraphMol()
                    mol.SetTitle(s)
                else:
                    mol.SetTitle(bio_number)
                    OESetSDData(mol,"BIO-NUMBER", bio_number)
            else:
                if mol is not None:
                    OESetSDData(mol,colNames[colId],str(row[colId].value))
        if mol is not None:
            OEWriteMolecule(ofs,mol)

    ofs.close()

