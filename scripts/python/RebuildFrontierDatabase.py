#!/usr/bin/env python

import psycopg2
import traceback
import os,sys
from socket import *
from openeye.oechem import *
from openeye.oemolprop import *
import requests
from bs4 import BeautifulSoup
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
        default_directory = "/DbUCdd/databases/VENDOR_LIBRARIES/"
    else:
        default_directory = "/Users/jfeng1/BiogenDB/Frontier/"


    try:
        login_url = "https://orders.frontierssi.com/FSSIOnline32/Login.aspx"
        download_url = "https://orders.frontierssi.com/FSSIOnline32/Downloads.aspx"
        session = requests.Session()
        headers={"User-Agent":"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_5) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/55.0.2883.95 Safari/537.36"}

        username = "jfeng"
        password = "temp123"

        #session.headers.update(headers)
        response1 = session.post(login_url)
        soup = BeautifulSoup(response1.content,"html.parser")

        VIEWSTATE=soup.find(id="__VIEWSTATE")['value']
        VIEWSTATEGENERATOR=soup.find(id="__VIEWSTATEGENERATOR")['value']
        EVENTVALIDATION=soup.find(id="__EVENTVALIDATION")['value']

        login_data={"__LASTFOCUS":"",
                    "__EVENTTARGET":"",
                    "__EVENTARGUMENT":"",
                    "__VIEWSTATE":VIEWSTATE,
                    "__VIEWSTATEGENERATOR":VIEWSTATEGENERATOR,
                    "__EVENTVALIDATION":EVENTVALIDATION,
                    "Login1$UserName":username,
                    "Login1$Password":password,
                    "Login1$LoginButton":"Log In"}

        session.post(login_url, data=login_data)

        response2 = session.post(download_url,verify=False)

        soup=BeautifulSoup(response2.content,"html.parser")
        VIEWSTATE=soup.find(id="__VIEWSTATE")['value']
        VIEWSTATEGENERATOR=soup.find(id="__VIEWSTATEGENERATOR")['value']
        EVENTVALIDATION=soup.find(id="__EVENTVALIDATION")['value']


        download_data={
            "__EVENTTARGET":"ctl00$ContentPlaceHolder1$RI_Download",
            "__EVENTARGUMENT":"",
            "__VIEWSTATE":VIEWSTATE,
            "__VIEWSTATEGENERATOR":VIEWSTATEGENERATOR,
            "__EVENTVALIDATION":EVENTVALIDATION
        }

        response3 = session.post(download_url,data=download_data)
        database = os.path.join(default_directory,"Frontier.sdf")
        open(database,"w").write(response3.content)
        # print response3.headers
        # input = raw_input("%s is downloaded, continue? (yes/no)"%database)
        # if input != "yes":
        #     print "Exiting..."
        #     sys.exit(1)

        smilesList = []
        ifs = oemolistream()
        ifs.open(database)
        mol = OEGraphMol()
        while OEReadMolecule(ifs,mol):
            smiles = OEMolToSmiles(mol).split(" ")[0]
            acdno = OEGetSDData(mol,"asdi_catno").strip()
            if len(acdno)>0:
                smilesList.append("%s %s"%(smiles,acdno))
        ifs.close()
        smilesBuffer = "\n".join(smilesList)
        dict = calculateCLogP(smilesBuffer)

        ifs = oemolistream()
        ifs.open(database)
        mol = OEGraphMol()
        
        conn = psycopg2.connect(database='frontier',user='medchem',host='javelin.biogen.com',password='medchem')
        cursor = conn.cursor()
        cursor.execute("drop table if exists frontier_chemicals;")
        cursor.execute("drop schema if exists rdk cascade;")
        cursor.execute("create table if not exists frontier_chemicals (chemical_id serial primary key not null, molfile text, acdno char(30), mw real, psa real, clogp real);")
        cursor.execute("create schema rdk;")
        n = 0
        while OEReadMolecule(ifs,mol):
            # if n == limit and limit > 0:
            #     break
            n += 1
            sys.stderr.write("\r%d molecules processed."%n)
            sys.stderr.flush()
            acdno = OEGetSDData(mol,"asdi_catno").strip()
            on_hand = OEGetSDData(mol,"on_hand").strip()
            if on_hand == "0":
                continue
            if len(acdno)==0:
                acdno = "frontier_noname_%d"%n
            if mol.NumAtoms() > 0:
                mol.SetTitle(acdno)
                molfile = convertMolToString(mol)
                OEDeleteEverythingExceptTheFirstLargestComponent(mol)
                mw = OECalculateMolecularWeight(mol)
                psa = OEGet2dPSA(mol)
                clogp = -999
                if dict.has_key(acdno):
                    clogp = dict[acdno]
                cursor.execute("insert into frontier_chemicals (molfile,acdno,psa,mw,clogp) values (%s,%s,%s,%s,%s)",(molfile,acdno,psa,mw,clogp))

        sys.stderr.write("\n")
        cursor.execute("select * into rdk.mols from (select chemical_id,mol_from_ctab(molfile::cstring) m from frontier_chemicals) tmp where m is not null;")
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
