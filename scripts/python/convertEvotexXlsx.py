#!/usr/bin/env python
from openpyxl import load_workbook
import glob
import os,sys
from openeye.oechem import *
DIRECTORY = "/cluster/database/databases/Evotec/Registered/"
BIIB_SDF = "database.sdf"

if __name__ == "__main__":
    dict = {}
    ifs = oemolistream()
    ifs.open(BIIB_SDF)
    mol = OEGraphMol()
    while OEReadMolecule(ifs,mol):
        dict[str(mol.GetTitle())] = OEGraphMol(mol)
    ifs.close()
    print len(dict)

    ofs = oemolostream()
    ofs.SetFormat(OEFormat_SDF)
    ofs.open("output.sdf")
    files = glob.glob(os.path.join(DIRECTORY,"*.xlsx"))
    for f in files:
        print f
        wb = load_workbook(f)
        active_sheet = wb.get_active_sheet()
        colNames = []
        first_row = active_sheet.rows[0]
        for c in first_row:
            colName = c.value
            if colName is not None:
                colNames.append(str(colName))
        for row in active_sheet.rows[1:]:
            compound_id = row[0].value
            if compound_id is not None:
                compound_id = str(compound_id)
                if dict.has_key(compound_id):
                    mol = dict[compound_id]
                    for id,colName in enumerate(colNames):
                        value = row[id].value
                        if value is not None:
                            OESetSDData(mol,str(colName), str(value))
                        else:
                            OESetSDData(mol,str(colName),"")
                    OEWriteMolecule(ofs,mol)
                else:
                    print compound_id
    ofs.close()


