import xml.etree.ElementTree as ET
import requests


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
    print xml
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
# payload = searchByTargetTemplate%(caller_name,"*wnk*")

structure = """test
  SciTegic04131517322D

 28 32  0  0  0  0            999 V2000
    5.0011    7.5280    0.0000 N   0  0
    5.0000    6.7192    0.0000 C   0  0
    5.6988    6.3156    0.0000 N   0  0
    5.6970    7.9316    0.0000 C   0  0
    6.3964    7.5316    0.0000 C   0  0
    6.3997    6.7217    0.0000 C   0  0
    7.7981    6.7274    0.0000 N   0  0
    7.7948    7.5373    0.0000 C   0  0
    7.0917    7.9439    0.0000 C   0  0
    7.0872    8.7504    0.0000 C   0  0
    7.7821    9.1546    0.0000 C   0  0
    7.7780    9.9604    0.0000 C   0  0
    7.0768   10.3606    0.0000 C   0  0
    6.3783    9.9492    0.0000 C   0  0
    6.3858    9.1448    0.0000 C   0  0
    7.0713   11.1671    0.0000 O   0  0
    7.7643   11.5769    0.0000 C   0  0
    8.4688   11.1768    0.0000 C   0  0
    9.1613   11.5860    0.0000 C   0  0
    9.1533   12.3919    0.0000 C   0  0
    8.4470   12.7869    0.0000 C   0  0
    7.7574   12.3755    0.0000 C   0  0
    8.3671    6.1593    0.0000 C   0  0
    5.6909    8.7419    0.0000 N   0  0
    8.2465    5.3644    0.0000 C   0  0
    8.9632    5.0000    0.0000 C   0  0
    9.5313    5.5691    0.0000 C   0  0
    9.1656    6.2851    0.0000 C   0  0
  1  2  2  0
  2  3  1  0
  3  6  2  0
  5  4  2  0
  4  1  1  0
  5  6  1  0
  7  8  1  0
  8  9  2  0
  9  5  1  0
  6  7  1  0
  9 10  1  0
 10 11  2  0
 11 12  1  0
 12 13  2  0
 13 14  1  0
 14 15  2  0
 15 10  1  0
 13 16  1  0
 16 17  1  0
 17 18  2  0
 18 19  1  0
 19 20  2  0
 20 21  1  0
 21 22  2  0
 22 17  1  0
  7 23  1  0
  4 24  1  0
 23 25  1  0
 25 26  1  0
 26 27  1  0
 27 28  1  0
 28 23  1  0
M  END
"""
searchByExactStructure = """<?xml version="1.0" encoding="UTF-8"?>
          <!DOCTYPE xf SYSTEM "rx.dtd">
          <xf>
            <request caller="%s" sessionid="">
              <statement command="select"/>
              <select_list>
                <select_item/>
              </select_list>
              <from_clause context="DPI" dbname="RX">
              </from_clause>
              <where_clause>structure('%s','compound,exact,isotopes,stereo_absolute,salts,mixtures,charges,radicals')</where_clause>
              <order_by_clause>DAT.PAUREUS desc</order_by_clause>
              <options>WORKER,NO_CORESULT</options>
            </request>
          </xf>"""

searchBySimilarity = """<?xml version="1.0" encoding="UTF-8"?>
          <!DOCTYPE xf SYSTEM "rx.dtd">
          <xf>
            <request caller="%s" sessionid="">
              <statement command="select"/>
              <select_list>
                <select_item/>
              </select_list>
              <from_clause context="S" dbname="RX">
              </from_clause>
              <where_clause>structure('%s','compound,similarity=60')</where_clause>
              <options>WORKER,NO_CORESULT</options>
            </request>
          </xf>"""

payload = searchByExactStructure%(caller_name,structure)

# And here you take the hitset resultname and you retrieve com.biogen.application.insilicotools.data from it, e.g. DAT (bioactivities), RX (reactions) etc.

# select_item_template = """                <select_item>%s</select_item>\n"""
# select_template = """<?xml version="1.0" encoding="UTF-8"?>
#           <!DOCTYPE xf SYSTEM "rx.dtd">
#           <xf>
#             <request caller="%s" sessionid="">
#               <statement command="select"/>
#               <select_list>\n"""
# for index in range (0,len(select_items)):
#     select_template = select_template + select_item_template%(select_items[index])
# select_template = select_template + """              </select_list>
#               <from_clause resultname="%s" %s first_item="%s" last_item="%s">
#               </from_clause>
#               <order_by_clause>%s</order_by_clause>
#               <group_by_clause>%s</group_by_clause>
#               <options>%s</options>
#             </request>
#           </xf>\n"""
# payload = select_template%(callername, resultname, grouplist, first_item, last_item, order_by, group_by, options)

# xml = login_response.text
# dict = getSessionInfo(xml)

resultXmlTemplate = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE xf SYSTEM "rx.dtd">
<xf>
<request caller="biogen_jun2015" sessionid="">
<statement command="select"/>
<select_list>
<select_item>IDE</select_item>
<select_item>DAT</select_item>
</select_list>
<from_clause resultname="%s"  first_item="1" last_item="10">
</from_clause>
<order_by_clause></order_by_clause>
<group_by_clause></group_by_clause>
<options></options>
</request>
</xf>"""
from openeye.oechem import *
ifs = oemolistream()
ifs.open("query.sdf")
mol = OEGraphMol()
molUtilities = MolUtilities()

ofs = oemolostream()
ofs.open("reaxys_output.sdf")
ofs.SetFormat(OEFormat_SDF)
while OEReadMolecule(ifs,mol):
    mdl = molUtilities.molToMDL(mol)
    payload = searchByExactStructure%(caller_name,mdl)
    print payload
    r = requests.post(reaxys_url,payload,cookies=cookies)
    xml = r.text.encode("utf-8")
    print xml
    resultNameAndSize = getResultNameAndSize(xml)
    print mol.GetTitle(),
    if resultNameAndSize is not None:
        print resultNameAndSize['result_size']
        payload = resultXmlTemplate%(resultNameAndSize['result_name'])
        r = requests.post(reaxys_url, payload, cookies=cookies)
        print r.text
        bioacts = getBioactivities(r.text.encode("utf-8"))
        for a in bioacts:
            mol1 = OEGraphMol(mol)
            a.attachToMol(mol1)
            OEWriteMolecule(ofs,mol1)
    else:
        print 0

ofs.close()
ifs.close()






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
