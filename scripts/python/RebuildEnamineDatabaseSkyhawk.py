#!/usr/bin/env python

import psycopg2
import traceback
from socket import *
from openeye.oechem import *
import requests
import zipfile
import os
import shutil
import xmlrpc
import json
from rdkit.Chem import Descriptors
from rdkit import Chem
import sys
from psycopg2 import sql

def calculateCLogP_PSA(smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        return logp,tpsa
    except:
        return 0,0

def convertMolToString(mol):
    ofs = oemolostream()
    ofs.SetFormat(OEFormat_SDF)
    ofs.openstring()
    OEWriteMolecule(ofs,mol)
    return ofs.GetString()

# limit = -1
if __name__ == "__main__":
    default_directory = None
    if gethostname() == "compchem":
        default_directory = "/home/jfeng/Downloads/"
    else:
        default_directory = "/Users/junfeng/Downloads/"

    #url = "http://www.enamine.net/index.php?option=com_content&task=view&id=10&auth=0"
    #url = "https://enamine.net/component/download/?f=4"
    #session = requests.Session()
    #response = session.post(url,data={"login":"jun@skyhawktx.com","password":"Beard301"})

    #session = requests.Session()
    #response = session.get("http://www.enamine.net/files/Building_Blocks/Enamine_BBstock_all.zip")
    #tmpfile = os.path.join(default_directory,"tmp.zip")
    #open(tmpfile,"wb").write(response.content)
    tmpfile = os.path.join(default_directory,"Enamine_BB.zip")
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
    #database = os.path.join(default_directory,"enamine_test.sdf")
#    input = raw_input("%s is downloaded, continue? (yes/no)"%database)
#    if input != "yes":
#        print "Exiting..."
#        sys.exit(1)

    try:
        smilesList = []
        ifs = oemolistream()
        ifs.open(database)
        logp_dict = {}
        psa_dict = {}
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            smiles = OEMolToSmiles(mol).split(" ")[0]
            logp,psa = calculateCLogP_PSA(smiles)
            catalogue_id = OEGetSDData(mol,"ID").strip()
            if len(catalogue_id)>0:
                smilesList.append("%s %s"%(smiles,catalogue_id))
                logp_dict[catalogue_id] = logp
                psa_dict[catalogue_id] = psa
        ifs.close()

        ifs = oemolistream()
        ifs.open(database)
        mol = OEGraphMol()
        
        conn = psycopg2.connect(database='Enamine',user='medchem',host='10.74.2.128',password='medchem')
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
            catalogue_id = OEGetSDData(mol,"ID").strip() #use "id" for catalogue sdf
            if len(catalogue_id)==0:
                continue
            if mol.NumAtoms() > 0:
                mol.SetTitle(catalogue_id)
                molfile = convertMolToString(mol).decode('UTF-8')
                OEDeleteEverythingExceptTheFirstLargestComponent(mol)
                mw = OECalculateMolecularWeight(mol)
                psa = 0.0
                if catalogue_id in psa_dict:
                    psa = psa_dict[catalogue_id]
                clogp = 0.0
                if catalogue_id in logp_dict:
                    clogp = logp_dict[catalogue_id]
                cursor.execute(sql.SQL("insert into enamine_chemicals (molfile,catalogue_id,psa,mw,clogp) values (%s,%s,%s,%s,%s)"),[molfile,catalogue_id,psa,mw,clogp])
                #cursor.execute("insert into enamine_chemicals (molfile,catalogue_id,psa,mw,clogp) values (%s,%s,%s,%s,%s)",(molfile,catalogue_id,psa,mw,clogp))
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
