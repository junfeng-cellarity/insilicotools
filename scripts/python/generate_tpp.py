#!/usr/bin/env python

import requests
import json
import time
import psycopg2
import numpy as np
from openeye.oechem import *
import progressbar
import sys


def molFromString(mol_string):
    oemol = OEGraphMol()
    ifs = oemolistream()
    ifs.SetFormat(OEFormat_SDF)
    ifs.openstring(mol_string)
    OEReadMolecule(ifs,oemol)
    ifs.close()
    return oemol

class SkyhawkDb:
    def __init__(self):
        self.conn = psycopg2.connect(database='skyhawk',user='medchem',host='10.74.2.128',password='medchem')
        
    def get_molecule(self,id):
        cursor = None
        try:
            cursor = self.conn.cursor()
            cursor.execute("select molfile,corporate_id from chemicals where chemical_id = %s",(id,))
            row = cursor.fetchone()
            if row is not None:
                return row[0],row[1]
            return None
        finally:
            if cursor is not None:
                cursor.close()
                
    def close(self):
        self.conn.close()
        
    def __del__(self):
        self.conn.close()
        

class ReadoutDef:
    def __init__(self, jsonText):
        dict = json.loads(jsonText)
        self.name = dict['name']
        self.data_type = dict['data_type']
        self.unit_label = dict['nM']
        self.id = dict['id']
    def __str__(self):
        return "%d %s %s %s"%(self.id,self.name,self.data_type,self.unit_label)

adme_protocol_names = [
    "General Metabolic Stability",
    "General hERG",
    "General MDCK-MDR1",
    "General Solubility",
    "General MDCK-MDR1 (Pharmaron Experimental Protocol)",
    "General Nova Ion Channels",
    "All Full PK",
    "All Neuro PK",
    "General Protein Binding"
]
       
keywords = ["data_sets","molecules","projects","protocols"]
api_dict = json.loads(requests.request("GET","http://10.74.2.128:8080/data/cdd_api_key.json").text)
cdd_url = "https://app.collaborativedrug.com/api/v1/vaults/4647"
headers = {'X-CDD-token':api_dict['production']}
mol_url = "https://app.collaborativedrug.com/vaults/4647/molecules/%s"

def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

def is_async_export_finished(id):
    dict = json.loads(requests.request("GET","%s/export_progress/%d"%(cdd_url,id),headers=headers).text)
    if dict['status']=="finished":
        return True
    else:
        return False

def get_cdd_protocols():
    dict = json.loads(requests.request("GET","%s/protocols/?page_size=1000"%(cdd_url), headers=headers).text)
    return dict['objects']

def get_protocol_data(id):
    max_time = 100
    dict = json.loads(requests.request("GET","%s/protocols/%d/data?async=true"%(cdd_url,id),headers=headers).text)
    export_id = dict['id']
    success = False
    start_time = time.time()
    pgbar = progressbar.ProgressBar(0,100)
    while True:
        time.sleep(1)
        elapsed_time = time.time()-start_time+1
        pgbar.update(int(elapsed_time))
        if elapsed_time>max_time:
            break
        if is_async_export_finished(export_id):
            success = True
            break
    if success:
        return (requests.request("GET","%s/exports/%d"%(cdd_url,export_id),headers=headers)).text
    else:
        return None


protocols = get_cdd_protocols()
project_protocols = []
adme_protocols = []


project_dict = {
    "sma":["SMN2 Potency"]
}

project_name = "sma"

protocol_prefixes = project_dict[project_name]
protocol_names = []
for object in protocols:
    for protocol_prefix in protocol_prefixes:
        if object['name'].startswith(protocol_prefix):
            project_protocols.append(object)
            protocol_names.append(object['name'])
            break
    if object['name'] in adme_protocol_names:
        adme_protocols.append(object)

project_data_dict = {}
project_definitions_dict = {}
for project_protocol in project_protocols:
    protocol_name = project_protocol['name']
    tmp_data = get_protocol_data(project_protocol['id'])
    project_data_dict[protocol_name] = tmp_data
    project_definitions_dict[protocol_name] = project_protocol['readout_definitions']

adme_data_dict = {}
adme_definitions_dict = {}
for adme_protocol in adme_protocols:
    protocol_name = adme_protocol['name']
    tmp_data = get_protocol_data(adme_protocol['id'])
    adme_data_dict[protocol_name] = tmp_data
    adme_definitions_dict[protocol_name] = adme_protocol['readout_definitions']
#print(adme_definitions_dict)
#print(adme_data_dict)

#{'molecule': 53706260, 'readouts': {'413437': {'modifier': '>', 'note': '(n=4)', 'value': 10000.0}, '414205': {'modifier': '<=>', 'value': 1.0}, '413440': {'modifier': '<=>', 'value': 1.0}, '412001': {'modifier': '<=>', 'value': 1.0}, '412003': {'modifier': '>', 'value': 3784.4}, '411999': {'modifier': '>', 'value': 10000.0}, '413438': {'modifier': '>', 'note': '(n=4)', 'value': 10000.0}, '412000': {'modifier': '>', 'value': 10000.0}}, 'batch': 51018110, 'id': 264405693, 'run': 235896}
def_dict = {}
for protocol_name in protocol_names:
    print (protocol_name)
    definitions = project_definitions_dict[protocol_name]
    for definition in definitions:
        def_dict[str(definition['id'])] = definition
for adme_name in adme_protocol_names:
    if len(adme_name)>0:
        definitions = adme_definitions_dict[adme_name]
        for definition in definitions:
            def_dict[str(definition['id'])] = definition
#print(def_dict)

mol_dict = {}
mol_id_list = []
i = 0
mol_count = 0
mol_data_dict = {}
mol_readout_dict = {}
for protocol_name in project_data_dict.keys():
    data_dict = json.loads(project_data_dict[protocol_name])
    molecules = data_dict['objects']
    for molecule in molecules:
        mol_id = molecule['molecule']
        if mol_id not in mol_id_list:
            mol_id_list.append(mol_id)
            mol_data_dict[mol_id] = {}


for mol_id in mol_id_list:
    for protocol_name in protocol_names:
        protocol_data = json.loads(project_data_dict[protocol_name])
        for mol in protocol_data['objects']:
            if mol['molecule'] == mol_id:
                for key in mol['readouts'].keys():
                    definition = def_dict[key]
                    def_name = definition['name']
                    def_type = definition['data_type']
                    if def_type == 'Number':
                        data_point = mol['readouts'][key]
                        if def_name not in mol_data_dict[mol_id]:
                            mol_data_dict[mol_id][def_name] = []
                        if isinstance(data_point,(int,float)):
                            mol_data_dict[mol_id][def_name].append(data_point)
                        elif 'value' in data_point:
                            mol_data_dict[mol_id][def_name].append(data_point['value'])
                        else:
                            continue


    #specis_id = '416380'
    pgp_inhibitor_id = '418460'
    for adme in adme_protocol_names:
        adme_data = json.loads(adme_data_dict[adme])
        for mol in adme_data['objects']:
            if mol['molecule'] == mol_id:
                specis = None
                specis_id = -1
                for key in mol['readouts'].keys():
                    if def_dict[key]['name'] == "Species":
                        specis_id = key
                        break
                if specis_id in mol['readouts'].keys():
                    specis = mol['readouts'][specis_id]

                RoA = None
                RoA_id = -1
                for key in mol['readouts'].keys():
                    if def_dict[key]['name']=="RoA":
                        RoA_id = key
                        break
                if RoA_id in mol['readouts'].keys():
                    RoA = mol['readouts'][RoA_id]

                pgp_inhibitor = None
                if pgp_inhibitor_id in mol['readouts'].keys():
                    pgp_inhibitor = mol['readouts'][pgp_inhibitor_id]
                for key in mol['readouts'].keys():
                    definition = def_dict[key]
                    def_name = "%s:%s"%(adme,definition['name'])
                    def_type = definition['data_type']
                    if def_type == 'Number':
                        if specis is not None:
                            def_name = "%s_%s" % (def_name, specis)
                        if pgp_inhibitor is not None and pgp_inhibitor == "YES":
                            def_name = "%s_%s" % (def_name, pgp_inhibitor)
                        if RoA is not None:
                            def_name = "%s_%s"%(def_name, RoA)
                        data_point = mol['readouts'][key]
                        if def_name not in mol_data_dict[mol_id]:
                            mol_data_dict[mol_id][def_name] = []
                        if isinstance(data_point,(int,float)):
                            mol_data_dict[mol_id][def_name].append(data_point)
                        elif 'value' in data_point:
                            mol_data_dict[mol_id][def_name].append(data_point['value'])
                        else:
                            continue
molDb = SkyhawkDb()
ofs = oemolostream()
ofs.open("/home/jfeng/apache-tomcat-9.0.13/webapps/data/%s_tpp2.sdf"%project_name)
desired_columns = ["Protein FL EC50","EC[1.5 Fold]", "CTG CC50", "Splicing FL EC50", "Splicing D7 IC50",
                   "FHOD3 IC50 Off-target", "Selectivity Marker (FHOD3/D7)", "SENP6 IC50 Off-target",
                   "Cytotoxicity marker (SENP6/D7)", "General MDCK-MDR1:Efflux", "General MDCK-MDR1 (Pharmaron Experimental Protocol):Efflux","General MDCK-MDR1:Papp A-B","General MDCK-MDR1 (Pharmaron Experimental Protocol):Papp A-B", "General hERG:hERG IC50", "General hERG:hERG IC50 Pharmaron",
                   "General Nova Ion Channels:hERG IC50",
                   "General Solubility:SGF (uM)", "General Solubility:FaSSIF (pH 6.5) (uM)",
                   "General Solubility:PBS (pH 7.4) (uM)", "General Metabolic Stability:T1/2_human","General Metabolic Stability:T1/2_mouse",
                   "Tier 1 Separator",
                   "General Protein Binding:Fu,b_mouse","General Protein Binding:Fu,p_mouse","General Protein Binding:Fu,p_human",
                   "All Full PK:Kp,uu AUC (BALB/c mouse)_Mouse",
                   "All Full PK:CL_Mouse_IV","All Full PK:CL_Mouse_PO","All Full PK:CL_Mouse_SC",
                   "All Full PK:AUCinf,b_Mouse_IV","All Full PK:AUCinf,b_Mouse_PO","All Full PK:AUCinf,b_Mouse_SC",
                   "All Full PK:AUCinf,p_Mouse_IV","All Full PK:AUCinf,p_Mouse_","All Full PK:AUCinf,p_Mouse_SC",
                   "Tier 2 Separator",
                   "General Nova Ion Channels:Cav1.2 IC50","General Nova Ion Channels:KvLQT1 IC50","General Nova Ion Channels:Late Nav1.5 IC50"
                   ]
smarts_dict = {
    "W": "F[#6]-1-[#6]-[#7]-[#6]-[#6]-[#6]-1-[#7]-c1nncs1",
    "Z": "F[#6]-1-[#6]-[#7]-[#6]-[#6]-[#6]-1-[#7]-c1cncnn1",
    "Y": "F[#6]-1-[#6]-[#7]-[#6]-[#6]-[#6]-1-[#7]-c1cnccn1"
}

for mol_id in mol_data_dict.keys():
    result = molDb.get_molecule(mol_id)
    if result is None:
        print("mol_id %d is not found, please sync the database."%mol_id)
        continue
    molfile,mol_name = result
    oemol = molFromString(molfile)
    oemol.SetTitle(mol_name)
    OESetSDData(oemol,"mol_id",str(mol_id))
    OESetSDData(oemol,"mol_url",mol_url%str(mol_id))
    found = False
    for series in smarts_dict.keys():
        ss = OESubSearch(smarts_dict[series])
        OEPrepareSearch(oemol, ss)
        if ss.SingleMatch(oemol):
            OESetSDData(oemol,"Series",series)
            found = True
            break
    if not found:
        OESetSDData(oemol,"Series","Other")

    # for tag in mol_data_dict[mol_id].keys():
    #         try:
    #             values = np.array(mol_data_dict[mol_id][tag])
    #             avg_value = np.mean(values)
    #             OESetSDData(oemol,tag,"%5.3f"%avg_value)
    #         except:
    #             OESetSDData(oemol,tag,"")

    for tag in desired_columns:
        if tag in mol_data_dict[mol_id].keys():
                try:
                    values = np.array(mol_data_dict[mol_id][tag])
                    avg_value = np.mean(values)
                    OESetSDData(oemol,tag,"%5.3f"%avg_value)
                except:
                    OESetSDData(oemol,tag,"")
        else:
            OESetSDData(oemol,tag,"")

    OEWriteMolecule(ofs,oemol)
ofs.close()



