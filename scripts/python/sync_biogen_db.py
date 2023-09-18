#!/usr/bin/env python
import cx_Oracle,base64
import os,datetime
from openeye.oechem import *
import psycopg2,traceback,sys

def convertMolToString(mol):
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_MDL)
    ofs.openstring()
    OEWriteMolecule(ofs,mol)
    return ofs.GetString()

def convertStringToMol(sdfString):
    ifs = oemolistream()
    ifs.SetFormat(OEFormat_SDF)
    ifs.openstring(sdfString)
    mol = OEGraphMol()
    OEReadMolecule(ifs,mol)
    ifs.close()
    return mol

#time = os.path.getmtime("/Users/jfeng1/biogen_addon.sdf")
#date = (datetime.datetime.fromtimestamp(time)-datetime.timedelta(days=1)).strftime("%d-%b-%y").upper()
#print date

timestamp_file = "/Users/jfeng1/.biogendb_timestamp"
date = open(timestamp_file,"r").read().strip()

date_now = datetime.datetime.now().strftime("%d-%b-%y").upper()
if date == date_now:
    print "Already updated today."
    sys.exit(0)

dbURL = base64.decodestring("amRiYzpvcmFjbGU6dGhpbjpALy8xMC4yLjEyOS4zNDoxNzI1L1BSRUxOUg==\n")
db_user = base64.decodestring("Y2hlbV91c2Vy\n")
db_password = base64.decodestring("QmlvZ2VuaWRlYzEyMw==\n")

connection = cx_Oracle.connect(db_user,db_password,"10.2.129.34:1725/PRELNR")
cursor = connection.cursor()
sql = "select compound_chemistry as MolValue, name as BIONumber from crdatamart.compounds where created_date > to_date('%s','DD-MON-YY')"%date
#sql = "select count(*) from crdatamart.compounds where created_date > to_date('%s','DD-MON-YY')"%date
cursor.execute(sql)
result = cursor.fetchall()
data = []
for molfile,bio_number in result:
    mol = convertStringToMol(str(molfile))
    mol.SetTitle(bio_number)                             
    data.append((mol, bio_number))
cursor.close()

print "%d new molecules found."%len(data)
if len(data)==0:
    sys.exit(0)
try:
    conn = psycopg2.connect(database='biogen2',user='medchem',host='javelin.biogen.com',password='medchem')
    cursor = conn.cursor()
    n = 0
    bios = {}

    for mol,bio_number in data:
        n += 1
        sys.stderr.write("\r%d molecules processed."%n)
        sys.stderr.flush()
        projname = "Unknown"
        if mol.NumAtoms() > 0 and not bios.has_key(bio_number):
            bios[bio_number] = 1
            mol.SetTitle(bio_number)
            molfile = convertMolToString(mol)
            cursor.execute("insert into chemicals (bio_number,molfile,project) select %s,%s,%s where not exists (select bio_number from chemicals where bio_number = %s)",(bio_number,molfile,projname,bio_number))
            cursor.execute("insert into rdk.mols (bio_number,m) select %s, mol_from_ctab(%s) where not exists (select bio_number from rdk.mols where bio_number = %s)",(bio_number,molfile,bio_number))
            cursor.execute("insert into rdk.fps (bio_number, torsionbv, mfp2, ffp2) select %s, torsionbv_fp(mol_from_ctab(%s)), morganbv_fp(mol_from_ctab(%s)), featmorganbv_fp(mol_from_ctab(%s)) where not exists (select bio_number from rdk.fps where bio_number = %s)",(bio_number,molfile,molfile,molfile,bio_number))
    sys.stderr.write("\n")
    cursor.close()
    conn.commit()

except:
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)

open(timestamp_file,"w").write(date_now)    
    

