#!/usr/bin/env python
from openpyxl import load_workbook
import glob
import os,sys
import psycopg2
import psycopg2.extras
import traceback
import numpy
from openeye.oechem import *
from openeye.oegraphsim import *
DIRECTORY = "/Users/jfeng1/Datasets/Fibroblast"

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
            traceback.print_exc()
            return None
        finally:
            if cursor is not None:
                cursor.close()

    def __del__(self):
        self.conn.close()


if __name__ == "__main__":
    db = BiogenDb()

    ifs = oemolistream()
    molDbPath = os.path.join(DIRECTORY,"BGNA8_proposed_plate_selection.sdf")
    ifs.open(molDbPath)
    molDb = OEMolDatabase(ifs)
    nrmols = molDb.GetMaxMolIdx()

    fpdb = OEFPDatabase(OEFPType_Path)

    emptyfp = OEFingerPrint()
    emptyfp.SetFPTypeBase(fpdb.GetFPTypeBase())

    for idx in range(0, nrmols):
        mol = OEGraphMol()
        if molDb.GetMolecule(mol, idx):
            fpdb.AddFP(mol)
        else:
            fpdb.AddFP(emptyfp)

    nrfps = fpdb.NumFingerPrints()

    plateDb = {}
    files = glob.glob(os.path.join(DIRECTORY,"*.xlsx"))
    for f in files:
        if os.path.basename(f).startswith("~"):
            continue
        wb = load_workbook(f)
        for sheet_name in wb.get_sheet_names():
            active_sheet = wb.get_sheet_by_name(sheet_name)
            progress = 0
            for rowid,row in enumerate(active_sheet.rows):
                plateName = None
                bio_number = None
                progress += 1
                sys.stderr.write("\r%d molecules processed."%progress)
                sys.stderr.flush()

                if sheet_name == "Result":
                    if rowid != 0:
                        plateName = row[0].value
                        bio_number = row[4].value

                else:
                    plateName = row[1].value
                    bio_number = row[7].value
                print progress
                if plateName is not None and bio_number is not None:
                    if not plateDb.has_key(plateName):
                        plateDb[plateName] = []
                    bio_number = bio_number.rsplit('-', 1)[0]
                    mol = db.getMolFromBio(bio_number)
                    if mol is not None:
                        hits = fpdb.GetSortedScores(mol, 1)
                        for hit in hits:
                            plateDb[plateName].append(hit.GetScore())

    plateNames = plateDb.keys()
    plates = []
    for plateName in plateNames:
        plates.append((plateName, numpy.average(plateDb[plateName])))
    sorted(plates,key=lambda plate:plate[1])
    f = open("plate_result.txt","w")
    for plateName,average in plates:
        print >> f, plateName,average
    f.close()




