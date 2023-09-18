import xml.etree.ElementTree as ET
import requests
from openeye.oechem import *

reaxys_url = "https://www.reaxys.com/reaxys/api"
caller_name = "biogen_jun2015"
DOI = "http://dx.doi.org/"
class MolUtilities():
    def __init__(self):
        return

    def convertToMol(self, molString, format):
        ifs = oemolistream()
        ifs.SetFormat(format)
        ifs.openstring(molString)
        mol = OEGraphMol()
        OEReadMolecule(ifs,mol)
        return mol

    def convertToInchiKey(self,mol):
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_INCHIKEY)
        ofs.openstring()
        OEWriteMolecule(ofs,mol)
        return ofs.GetString()

    def molToMDL(self, mol):
        ofs = oemolostream()
        ofs.SetFormat(OEFormat_MDL)
        ofs.openstring()
        OEWriteMolecule(ofs,mol)
        return ofs.GetString()

class Bioactivity:
    def __init__(self, refType, reference, targetName, vtype, value, unit, dose):
        self.refType = refType
        self.reference = reference
        if refType == "Article":
            self.reference = "%s%s"%(DOI,reference)
        self.targetName = targetName
        self.vtype = vtype
        try:
            self.value = float(value)
        except:
            self.value = None
        self.unit = unit
        self.dose = dose
        return

    def attachToMol(self,mol):
        if self.reference is not None:
            OESetSDData(mol, "Reference", self.reference)
        if self.targetName is not None:
            OESetSDData(mol, "Target", self.targetName)
        if self.vtype is not None:
            OESetSDData(mol, "Value Type", self.vtype)
        if self.value is not None:
            OESetSDData(mol, "Value", "%5.2f"%self.value)
        if self.unit is not None:
            OESetSDData(mol, "Unit", self.unit)
        if self.vtype == "% Inhibition":
            if self.dose is not None:
                OESetSDData(mol, "Dose", self.dose)

    def __str__(self):
        if self.vtype == "% Inhibition":
            return "%s %s %5.2f %s @%s %s"%(self.targetName,self.vtype,self.value,self.unit, self.dose, self.reference)
        else:
            return "%s %s %5.2f %s %s"%(self.targetName,self.vtype,self.value,self.unit,self.reference)

def getBioactivities(xml):
    bioactivities = []
    root = ET.fromstring(xml)
    for r in root.iter("dpitem"):
        refType = None
        reference = None
        targetName = None
        vtype = None
        value = None
        unit = None
        dose = None
        for v in r.iter("DAT.TNAME"):
            targetName = v.text
            break
        for v in  r.iter("DAT.VTYPE"):
            vtype = v.text
            break
        for v in r.iter("DAT.VALUE"):
            value = v.text
            break
        for u in r.iter("DAT.UNIT"):
            unit = u.text.encode("utf-8")
            break
        for u in r.iter("DAT.MDOSE"):
            dose = u.text.encode("utf-8")
            break

        for v in r.iter("CIT.DT"):
            refType = v.text

        if refType is not None:
            if refType == "Article":
                for v in r.iter("CIT.DOI"):
                    reference = v.text
                    break
            else:
                ptn = ""
                page = ""
                for v in r.iter("CIT.PN"):
                    ptn = v.text
                    break
                for v in r.iter("CIT.PK"):
                    page = v.text
                    break
                reference = "%s %s"%(ptn,page)
        bioactivity = Bioactivity(refType, reference, targetName, vtype, value, unit, dose)
        if bioactivity.value is not None and bioactivity.targetName is not None and bioactivity.reference is not None and bioactivity.unit is not None:
            bioactivities.append(bioactivity)
    return bioactivities

def getResultNameAndSize(xml):
    rootNode = ET.fromstring(xml)
    result_name = None
    result_size = 0
    for childNode in rootNode:
        if childNode.tag == "response":
            statusNode = childNode.find("status")
            if statusNode is None or statusNode.text != "OK":
                break
            results = childNode.find("results")
            if results is not None:
                resultNode = results.find("result")
                if resultNode is not None:
                    resultnameNode = resultNode.find("resultname")
                    if resultnameNode is not None:
                        result_name = resultnameNode.text
                    resultsizeNode = resultNode.find("resultsize")
                    if resultsizeNode is not None:
                        result_size = int(resultsizeNode.text)
    if result_name is None or result_size == 0:
        return None
    else:
        return {"result_name":result_name,"result_size":result_size}

def getSessionInfo(xml):
    rootNode = ET.fromstring(xml)
    stationid = None
    sessionid = None
    sessiontoken = None
    for childNode in rootNode:
        if childNode.tag == "response":
            statusNode = childNode.find("status")
            if statusNode is None or statusNode.text != "OK":
                break

            sessions = childNode.find('sessions')
            if sessions is None or len(sessions) == 0:
                break

            session = sessions[0]

            sessionidNode = session.find("sessionid")
            if sessionidNode is not None:
                sessionid = sessionidNode.text

            tokenNode = session.find("session_token")
            if tokenNode is not None:
                sessiontoken = tokenNode.text

            stationNode = session.find("stationid")
            if statusNode is not None:
                stationid = stationNode.text

    if stationid is not None and sessionid is not None and sessiontoken is not None:
        return {'stationid':stationid,'sessionid':sessionid,'sessiontoken':sessiontoken}

    return None


login_request = """<?xml version="1.0"?><!DOCTYPE xf SYSTEM "https://www.reaxys.com/xfserv/rx.dtd"><xf><request caller="%s"><statement command="connect" stationid="" ip_address="" username="" password=""/></request></xf>"""%caller_name
login_response=requests.post(reaxys_url,login_request)
cookies = login_response.cookies


# Description	Description of facts (XX) and fields (XX.XX)
# Code	Codes of facts (XX) and fields (XX.XX)
# Context	"Context to which a field (XX.XX) is associated with. Possible contexts are:
# Substances - S
# Reactions - R
# Citations - C
# Bioactivity Data Points - DPI
# Targets - TPI"
# Searchtype	"Description of whether and how a field (XX.XX) can be searched. Possible searchtypes are:
# S - searchable by text, e.g. IDE.CN='aspirin'
# I - searchable by number, e.g. RXD.T=20
# R - searchable by range, e.g. RXD.T=15-25
# A blank indicates that the field cannot be searched."
# Available for	"Description of the license that is required to search a field (XX.XX) and retrieve com.biogen.application.insilicotools.data from the field. Possible licenses are:
#     RX - Reaxys only license
# RMC - Reaxys Medicinal Chemistry only license
# combined - both Reaxys and Reaxys Medicinal Chemistry license
# RMC FF - field content is available from the Reaxys Medicinal Chemistry Flatfile"


searchByTargetTemplate = """<?xml version="1.0" encoding="UTF-8"?>
          <!DOCTYPE xf SYSTEM "rx.dtd">
          <xf>
            <request caller="%s" sessionid="">
              <statement command="select"/>
              <select_list>
                <select_item/>
              </select_list>
              <from_clause context="DPI" dbname="RX">
              </from_clause>
              <where_clause>TARGET.NAME='%s'</where_clause>
              <order_by_clause>DAT.PAUREUS desc</order_by_clause>
              <options>WORKER,NO_CORESULT</options>
            </request>
          </xf>"""
payload = searchByTargetTemplate%(caller_name,"PINK1")
r = requests.post(reaxys_url,payload,cookies=cookies)
print r.text.encode('utf-8').strip()
resultNameAndSize = getResultNameAndSize(r.text)

resultXmlTemplate = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xf SYSTEM "rx.dtd">
<xf>
<request caller="biogen_jun2015" sessionid="">
<statement command="select"/>
<select_list>
<select_item>DAT</select_item>
<select_item>IDE</select_item>
</select_list>
<from_clause resultname="%s"  first_item="1" last_item="%d">
</from_clause>
<order_by_clause></order_by_clause>
<group_by_clause></group_by_clause>
<options></options>
</request>
</xf>"""

payload = resultXmlTemplate%(resultNameAndSize['result_name'],resultNameAndSize['result_size'])
r = requests.post(reaxys_url, payload, cookies=cookies)

print r.text.encode('utf-8').strip()




# print dict
#
# # search_request = """<?xml version="1.0"?><!DOCTYPE xf SYSTEM "https://www.reaxys.com/xfserv/rx.dtd"><xf><request caller="biogen_jun2015" sessionid="%s" stationid="%s"><statement command="select"/><select_list><select_item></select_item><select_item>IDE</select_item>MP(1,10)</select_list><from_clause resultname="H001_123" first_item="20" last_item="20"/></request></xf>"""%(dict['sessionid'],dict['stationid'])
# #
# # req = urllib2.Request(reaxys_url,search_request)
# # search_response = urllib2.urlopen(req)
# # print search_response.read()
#


logout_request = """<?xml version="1.0"?><!DOCTYPE xf SYSTEM "https://www.reaxys.com/xfserv/rx.dtd"><xf><request caller="biogen_jun2015"><statement command="disconnect"/></request></xf>"""
r = requests.post(reaxys_url,logout_request,cookies=cookies)
