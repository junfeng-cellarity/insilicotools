#!/usr/bin/env python

import psycopg2
import traceback
import sys
from socket import *
from openeye.oechem import *
from openeye.oemolprop import *
import requests
import zipfile
import os
import shutil
import xmlrpclib
import json

def calculateCLogP(smilesBuffer):
    s = xmlrpclib.Server("http://javelin.corp.biogen.com:9528/")
    result = s.CLogPBatch(smilesBuffer)
    return json.loads(result)

def convertMolToString(mol):
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_MDL)
    ofs.openstring()
    OEWriteMolecule(ofs,mol)
    return ofs.GetString()

# limit = -1
if __name__ == "__main__":
    default_directory = None
    if gethostname() == "javelin" or gethostname() == "modeling":
        default_directory = "/DbUCdd/databases/BUILDING_BLOCKS/"
    else:
        default_directory = "/Users/jfeng1/BiogenDB/Enamine/"

    # url = "http://www.enamine.net/index.php?option=com_content&task=view&id=10&auth=0"
    # session = requests.Session()
    # response = session.post(url,data={"login":"jfeng1","password":"Beard301"})
    session = requests.Session()
    response = session.get("http://www.enamine.net/files/Building_Blocks/Enamine_BBstock_all.zip")
    tmpfile = os.path.join(default_directory,"tmp.zip")
    open(tmpfile,"wb").write(response.content)
    zip_file = zipfile.ZipFile(tmpfile)
    sdf_filename = None
    for name in zip_file.namelist():
        if name.endswith(".sdf"):
            sdf_filename = name
            break
    zip_file.extract(sdf_filename,path=default_directory)
    database = os.path.join(default_directory,sdf_filename)
    database2 = os.path.join(default_directory,"Enamine_BB_current.sdf")
    shutil.copyfile(database,database2)
#    input = raw_input("%s is downloaded, continue? (yes/no)"%database)
#    if input != "yes":
#        print "Exiting..."
#        sys.exit(1)

    try:
        smilesList = []
        ifs = oemolistream()
        ifs.open(database)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            smiles = OEMolToSmiles(mol).split(" ")[0]
            catalogue_id = OEGetSDData(mol,"id").strip()
            if len(catalogue_id)>0:
                smilesList.append("%s %s"%(smiles,catalogue_id))
        ifs.close()
        smilesBuffer = "\n".join(smilesList)
        dict = calculateCLogP(smilesBuffer)


        ifs = oemolistream()
        ifs.open(database)
        mol = OEGraphMol()
        
        conn = psycopg2.connect(database='Enamine',user='medchem',host='javelin.biogen.com',password='medchem')
        cursor = conn.cursor()
        cursor.execute("drop table if exists enamine_chemicals;")
        cursor.execute("drop schema if exists rdk cascade;")
        cursor.execute("create table if not exists enamine_chemicals (chemical_id serial primary key not null, molfile text, catalogue_id char(30), mw real, psa real, clogp real);")
        cursor.execute("create schema rdk;")
        n = 0
        while OEReadMolecule(ifs,mol):
            # if n == limit and limit > 0:
            #     break
            n += 1
            sys.stderr.write("\r%d molecules processed."%n)
            sys.stderr.flush()
            catalogue_id = OEGetSDData(mol,"IDNumber").strip() #use "id" for catalogue sdf
            if len(catalogue_id)==0:
                continue
            if mol.NumAtoms() > 0:
                mol.SetTitle(catalogue_id)
                molfile = convertMolToString(mol)
                OEDeleteEverythingExceptTheFirstLargestComponent(mol)
                mw = OECalculateMolecularWeight(mol)
                psa = OEGet2dPSA(mol)
                clogp = dict[catalogue_id]
                cursor.execute("insert into enamine_chemicals (molfile,catalogue_id,psa,mw,clogp) values (%s,%s,%s,%s,%s)",(molfile,catalogue_id,psa,mw,clogp))
        sys.stderr.write("\n")
        cursor.execute("select * into rdk.mols from (select chemical_id,mol_from_ctab(molfile::cstring) m from enamine_chemicals) tmp where m is not null;")
        cursor.execute("create index molidx on rdk.mols using gist(m);")
        cursor.execute("alter table rdk.mols add primary key (chemical_id);")
        cursor.execute("select chemical_id,torsionbv_fp(m) as torsionbv, morganbv_fp(m) as mfp2,featmorganbv_fp(m) as ffp2 into rdk.fps from rdk.mols;")
        cursor.execute("create index fps_ttbv_idx on rdk.fps using gist(torsionbv);")
        cursor.execute("create index fps_mfp2_idx on rdk.fps using gist(mfp2);")
        cursor.execute("create index fps_ffp2_idx on rdk.fps using gist(ffp2);")
        cursor.execute("alter table rdk.fps add primary key (chemical_id);")
        cursor.close()
        conn.commit()

    except:
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
