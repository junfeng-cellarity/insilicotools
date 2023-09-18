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


volumnColName = "Volume UL"
grossColName = "Gross Weight"
tareColName = "Tare Weight"
bionumberColName = "Corporate ID"
lotColName = "Lot Number"
barcodeName = "Barcode"


if __name__ == "__main__":
    if len(sys.argv)!=2:
        print "Usage:%s input.xls"%sys.argv[0]
        print "First row contains data tag names"
        sys.exit(0)

    OEThrow.SetLevel(OEErrorLevel_Error)
    xlsFile = sys.argv[1]

    liquidDb = {}
    solidDb = {}

    liquidBarcodeDb = {}
    solidBarcodeDb = {}

    wb = load_workbook(xlsFile)

    print "Loading done."
    solid_sheet = wb.get_sheet_by_name("Solid")
    colNames = []
    first_row = solid_sheet.rows[0]
    for c in first_row:
        colName = c.value
        colNames.append(colName.encode('utf-8').strip())

    try:
        bio_number_col_id = colNames.index(bionumberColName)
        grossId = colNames.index(grossColName)
        tareId = colNames.index(tareColName)
        lotId = colNames.index(lotColName)
        barcodeId = colNames.index(barcodeName)
    except:
        print "No desired columns"
        sys.exit(1)

    for row in solid_sheet.rows[1:]:
        try:
            bio_number = str(row[bio_number_col_id].value)
            lot_number = str(row[lotId].value)
            gross_weight = float(row[grossId].value)
            tare_weight = float(row[tareId].value)
            weight = 1000*(gross_weight-tare_weight)
            barcode = str(row[barcodeId].value)
            if weight <= 0:
                continue
            key = "%s-%s"%(bio_number,lot_number)
            if not solidDb.has_key(key):
                solidDb[key] = weight
                solidBarcodeDb[key] = []
                solidBarcodeDb[key].append(barcode)
            else:
                solidDb[key] = solidDb[bio_number]+weight
                solidBarcodeDb[key].append(barcode)
        except:
            traceback.print_exc()
            pass


    colNames = []
    liquid_sheet = wb.get_sheet_by_name("Liquid")
    first_row = liquid_sheet.rows[0]
    for c in first_row:
        colName = c.value
        colNames.append(colName.encode('utf-8').strip())

    bio_number_col_id = colNames.index(bionumberColName)
    volumColId = colNames.index(volumnColName)
    lotId = colNames.index(lotColName)
    barcodeId = colNames.index(barcodeName)
    for row in liquid_sheet.rows[1:]:
        try:
            bio_number = str(row[bio_number_col_id].value)
            lot_number = str(row[lotId].value)
            barcode = str(row[barcodeId].value)
            key = "%s-%s"%(bio_number,lot_number)
            volume = float(row[volumColId].value)
            if not liquidDb.has_key(key):
                liquidDb[key] = []
                liquidDb[key].append(volume)
                liquidBarcodeDb[key] = []
                liquidBarcodeDb[key].append(barcode)
            else:
                liquidDb[key].append(volume)
                liquidBarcodeDb[key].append(barcode)
        except:
            traceback.print_exc()
            pass

    n_keep = 0
    n_drop = 0
    drop_output = open("drop.txt","w")
    keep_output = open("keep.txt","w")
    for key in liquidDb.keys():
        if solidDb.has_key(key) and solidDb[key]>=2.0:
            if max(liquidDb[key])<150:
                n_drop += len(liquidDb[key])
                for id,barcode in enumerate(liquidBarcodeDb[key]):
                    print >> drop_output,barcode,key,liquidDb[key][id],"have_solid_%f"%solidDb[key]
            else:
                for id,v in enumerate(liquidDb[key]):
                    if v<150:
                        n_drop += 1
                        print >> drop_output,liquidBarcodeDb[key][id],key,v,"have_solid_%f"%solidDb[key]
                    else:
                        n_keep += 1
                        print >> keep_output,liquidBarcodeDb[key][id],key,v, "have_solid%f"%solidDb[key]
        else:
            if max(liquidDb[key])>=150:
                for id,v in enumerate(liquidDb[key]):
                    if v>=150:
                        n_keep += 1
                        print >> keep_output,liquidBarcodeDb[key][id],key,v,"liquid>=150_noSolid"
                    else:
                        n_drop += 1
                        print >> drop_output,liquidBarcodeDb[key][id],key,v,"liquid<150_noSolid"
            else:
                for id,v in enumerate(liquidDb[key]):
                    n_keep += 1
                    print >> keep_output,liquidBarcodeDb[key][id],key,v,"liquid<150_unique"


    drop_output.close()
    keep_output.close()
    print "No. liquid samples to keep:",n_keep
    print "No. liquid samples to drop:",n_drop
