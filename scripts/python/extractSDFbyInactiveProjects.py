from openpyxl import load_workbook
import glob
import os,sys
import psycopg2
import psycopg2.extras
import traceback
from parse import *
from openeye.oechem import *

class BiogenDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='biogen',user='medchem',host='javelin',password='medchem')

    def getMolFromBio(self, bio_number):
        cursor = None
        try:
            cursor = self.conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
            cursor.execute("select molfile,project from chemicals where bio_number = %s and project not in ('NRF2','IRAK4','BTK','SIK','FAK','PKR','GLK','RIPK1','MALT1','PLD1')",(bio_number,))
            result = cursor.fetchone()
            if result is not None:
                molfile = result['molfile']
                mol = OEGraphMol()
                ifs = oemolistream()
                ifs.SetFormat(OEFormat_MDL)
                ifs.openstring(molfile)
                OEReadMolecule(ifs,mol)
                mol.SetTitle(bio_number)
                OESetSDData(mol,"project",result['project'])
                return mol

            else:
                return None
        except:
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()


if __name__ == "__main__":
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_SDF)
    ofs.open("output.sdf")
    excel_file = "/Users/jfeng1/Biogen Inventory_9_21_15.xlsx"
    wb = load_workbook(excel_file)
    active_sheet = wb.get_active_sheet()
    colNames = []
    first_row = active_sheet.rows[0]
    db = BiogenDb()
    for c in first_row:
        colName = c.value
        if colName is not None:
            colNames.append(str(colName))
    for row in active_sheet.rows[1:]:
        compound_id = row[0].value
        onholdSolidAmount = row[1].value
        onholdLiquidAmount = row[2].value
        if compound_id is not None:
            compound_id = str(compound_id)
            new_id = None
            a = parse("BIO-{}-{}",compound_id)
            if a == None:
                a = parse("BIO-{}",compound_id)
                if a != None:
                    new_id = "BIO-%s"%(a[0])
            else:
                new_id = "BIO-%s"%(a[0])
            if new_id is not None:
                mol = db.getMolFromBio(new_id)
                if mol is not None:
                    OESetSDData(mol,colNames[0],compound_id)
                    OESetSDData(mol,colNames[1],str(onholdSolidAmount))
                    OESetSDData(mol,colNames[2],str(onholdLiquidAmount))
                    OEWriteMolecule(ofs,mol)
    ofs.close()


