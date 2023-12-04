__author__ = 'jfeng1'
import sys
import os
import traceback, tempfile, urllib2
from time import sleep
import javax.swing.ProgressMonitor as ProgressMonitor
import javax.swing.SwingWorker as SwingWorker
import java.net.URL as URL
import java.lang.Math as Math
import javax.swing.JComboBox as JComboBox
import javax.swing.JLabel as JLabel
import javax.swing.JPanel as JPanel
import java.util.Properties as Properties

import java.util.concurrent.ExecutionException as ExecutionException
import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.table.VortexJTable2 as VortexJTable2
import com.dotmatics.vortex.util.Util as Util
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
from com.dotmatics.chemlink import ChemLink
import com.dotmatics.vortex.table.IntArray as IntArray
import javax.swing.SwingUtilities as SwingUtilities
import java.lang.String as String
import org.python.core.util.StringUtil as StringUtil
import javax.swing.event.ChangeListener as ChangeListener
from javax.swing import BorderFactory, ProgressMonitor, SwingWorker, JScrollPane, JComboBox, JPanel, JTable, ImageIcon, \
    JLabel, JFileChooser, JCheckBox, JSlider, JToolBar
from javax.swing import JButton, JFrame, JTextField, BoxLayout
from javax.swing.border import TitledBorder
from java.awt import GridBagLayout as awtGridBagLayout, BorderLayout, FlowLayout, Dimension, GridLayout
from java.awt import GridBagConstraints as awtGridBagConstraints
from java.awt import Insets as awtInsets
from com.ziclix.python.sql import zxJDBC
from javax.swing.table import DefaultTableModel, DefaultTableCellRenderer
import jarray
import java.awt.image.BufferedImage as BufferedImage
import java.io.ByteArrayOutputStream as ByteArrayOutputStream
import java.awt.event.ActionListener as ActionListener
import java.awt.event.KeyListener as KeyListener
import java.awt.event.KeyEvent as KeyEvent
import javax.swing.JOptionPane as JOptionPane

import java.beans.PropertyChangeListener as PropertyChangeListener
sys.path.append(os.path.join(vortex.getVortexFolder(),"jars"))
from init_vortex import VortexInit
vortexInit = VortexInit(vortex)
vortexInit.initJars()

import chemaxon.marvin.beans.MSketchPane as MSketchPane
import chemaxon.formats.MolImporter as MolImporter
import chemaxon.formats.MolExporter as MolExporter

from chemaxon.formats import *
from chemaxon.struc import *
from chemaxon.marvin.plugin.CalculatorPlugin import *
from chemaxon.marvin.calculations import *
import chemaxon.marvin.calculations.pKaPlugin as PkaPlugin

import org.apache.xmlrpc.client.XmlRpcClient as XmlRpcClient
import org.apache.xmlrpc.client.XmlRpcClientConfigImpl as XmlRpcClientConfigImpl
import urllib, subprocess
import org.json.simple.JSONArray as JSONArray
import org.json.simple.JSONObject as JSONObject
import org.json.simple.parser.JSONParser as JSONParser
import org.json.simple.parser.ParseException as ParseException
import httplib
class MolUtilities:
    def __init__(self):
        return

    def convertToMol(self,molString):
        mol = MolImporter.importMol(molString,"mol")
        return mol

    def molToSmiles(self,mol):
        return MolExporter.exportToFormat(mol,"cxsmiles:u,a0")
    def convertMolsToSDF(self,mol_list):
        stream = ByteArrayOutputStream()
        molExporter = MolExporter(stream,"sdf")
        for id, mol in enumerate(mol_list):
            #mol.setProperty("moka_id","%s"%id)
            mol.setName("%d"%id)
            molExporter.write(mol)
        return stream.toString()
class ADMERunner:
    def __init__(self):
        self.columnTypeDict = {
            "ALogD":vortex.DOUBLE,
            "MoKa LogP":vortex.DOUBLE,
            "MoKa LogD":vortex.DOUBLE,
            "MoKa Basic pKa":vortex.DOUBLE,
            "MoKa Acidic pKa":vortex.DOUBLE,
            "ChemAxon LogP":vortex.DOUBLE,
            "ChemAxon LogD":vortex.DOUBLE,
            "ChemAxon Basic pKa":vortex.DOUBLE,
            "ChemAxon Acidic pKa":vortex.DOUBLE,
            "Cellarity_hERG_IC50(uM)":vortex.DOUBLE,
            "Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)":vortex.DOUBLE,
            "hepatocyte stability":vortex.STRING
        }
        self.properties = ["ALogD","MoKa LogP","MoKa LogD", "MoKa Basic pKa", "MoKa Acidic pKa",
                           "ChemAxon LogP", "ChemAxon LogD", "ChemAxon Basic pKa", "ChemAxon Acidic pKa",
                           "Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)",
                           "Cellarity_hERG_IC50(uM)", "hepatocyte stability"
                           ]
        self.xmlrpc_client = XmlRpcClient()
        self.xmlrpc_config = XmlRpcClientConfigImpl()
        self.xmlrpc_config.setServerURL(URL("http://10.74.2.128:9527"))
        self.xmlrpc_client.setConfig(self.xmlrpc_config)
        self.checkboxDict = {}

        self.pkaPlugin = pKaPlugin()
        pkaParams = Properties()
        pkaParams.put("type","PKa,acidic,basic")
        self.pkaPlugin.setParameters(pkaParams)
        self.pkaPlugin.setpH(7.4)
        self.pkaPlugin.setBasicpKaLowerLimit(-5)
        self.pkaPlugin.setAcidicpKaUpperLimit(20)
        self.logPlugin = logPPlugin()
        self.logDPlugin = logDPlugin()
        self.logDPlugin.setpH(7.4)

    def calculate_cLogD74_chemaxon(self,mol):
        plugin = self.logDPlugin
        plugin.setMolecule(mol)
        plugin.run()
        return plugin.getlogD()

    def calculate_cLogP_chemaxon(self,mol):
        plugin = self.logPlugin
        plugin.setMolecule(mol)
        plugin.run()
        return plugin.getlogPTrue()

    def calculate_pKa_chemaxon(self,mol):
        plugin = self.pkaPlugin
        plugin.setMolecule(mol)
        plugin.run()
        acidic_pka = plugin.getMacropKaValues(PkaPlugin.ACIDIC)
        basic_pka = plugin.getMacropKaValues(PkaPlugin.BASIC)
        return acidic_pka,basic_pka



    def build_property_panel(self,title,properties, ncolumn):
        nrow = len(properties)/2+1
        p = JPanel(GridLayout(nrow, ncolumn))
        p.setBorder(TitledBorder(title))
        for property in properties:
            cb = JCheckBox(property, True)
            p.add(cb)
            self.checkboxDict[property] = cb
        return p

    def calculate_hERG(self, smilesDict):
        conn = httplib.HTTPConnection("herg.cellarity.int", port=80)
        headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
        parser = JSONParser()

        smilesList = JSONArray()
        for key in smilesDict.keys():
            smilesEntry = JSONObject()
            smilesEntry.put("id", key)
            smilesEntry.put("smiles", smilesDict[key])
            smilesList.add(smilesEntry)
        conn.request("POST", "/predict", smilesList.toJSONString(), headers)
        response = conn.getresponse()
        if response.status == 200:
            data = response.read()
            result_dict = parser.parse(data)
            return result_dict
        return {}

    def calculate_Perm(self, smilesDict):
        conn = httplib.HTTPConnection("perm.cellarity.int", port=80)
        headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
        parser = JSONParser()

        smilesList = JSONArray()
        for key in smilesDict.keys():
            smilesEntry = JSONObject()
            smilesEntry.put("id", key)
            smilesEntry.put("smiles", smilesDict[key])
            smilesList.add(smilesEntry)
        conn.request("POST", "/predict", smilesList.toJSONString(), headers)
        response = conn.getresponse()
        if response.status == 200:
            data = response.read()
            result_dict = parser.parse(data)
            return result_dict
        else:
            print(response.msg)
            print("Something wrong.")
        return {}

    def calculate_MoKa(self, smilesDict):
        conn = httplib.HTTPConnection("moka.cellarity.int", port=80)
        headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
        parser = JSONParser()

        smilesList = JSONArray()
        for key in smilesDict.keys():
            smilesEntry = JSONObject()
            smilesEntry.put("id", key)
            smilesEntry.put("smiles", smilesDict[key])
            smilesList.add(smilesEntry)
        conn.request("POST", "/predict", smilesList.toJSONString(), headers)
        response = conn.getresponse()
        if response.status == 200:
            data = response.read()
            result_dict = parser.parse(data)
            return result_dict
        else:
            print(response.msg)
            print("Something wrong.")
        return {}

    def calculate_HepatocyteStability(self, smilesDict):
        conn = httplib.HTTPConnection("hepa.cellarity.int", port=80)
        headers = {'accept': 'application/json', 'Content-Type': 'application/json'}
        parser = JSONParser()

        smilesList = JSONArray()
        for key in smilesDict.keys():
            smilesEntry = JSONObject()
            smilesEntry.put("id", key)
            smilesEntry.put("smiles", smilesDict[key])
            smilesList.add(smilesEntry)
        conn.request("POST", "/predict", smilesList.toJSONString(), headers)
        response = conn.getresponse()
        if response.status == 200:
            data = response.read()
            result_dict = parser.parse(data)
            return result_dict
        else:
            print(response.msg)
            print("Something wrong.")
        return {}

    def show_interface(self):
        this = self
        mypanel = JPanel(BorderLayout())
        toolbar = JToolBar()
        selectAllBtn = JButton("Select All")
        class SelectListener(ActionListener):
            def actionPerformed(self,e):
                for key in this.checkboxDict.keys():
                    this.checkboxDict[key].setSelected(True)
        selectAllBtn.addActionListener(SelectListener())

        unselectAllBtn = JButton("Unselect All")
        class UnselectListener(ActionListener):
            def actionPerformed(self,e):
                for key in this.checkboxDict.keys():
                    this.checkboxDict[key].setSelected(False)
        unselectAllBtn.addActionListener(UnselectListener())

        invertSelectBtn = JButton("Invert Selection")
        class InvertListener(ActionListener):
            def actionPerformed(self,e):
                for key in this.checkboxDict.keys():
                    if this.checkboxDict[key].isSelected():
                       this.checkboxDict[key].setSelected(False)
                    else:
                        this.checkboxDict[key].setSelected(True)
        invertSelectBtn.addActionListener(InvertListener())

        toolbar.add(selectAllBtn)
        toolbar.add(unselectAllBtn)
        toolbar.add(invertSelectBtn)
        mypanel.add(toolbar,BorderLayout.NORTH)
        p1 = self.build_property_panel("Properties",self.properties,3)
        #p2 = self.build_property_panel("ADME",self.adme_models,2)
        mainPanel = JPanel()
        layout = BoxLayout(mainPanel,BoxLayout.Y_AXIS)
        mainPanel.setLayout(layout)
        mainPanel.add(p1)
        #mainPanel.add(p2)
        mypanel.add(JScrollPane(mainPanel),BorderLayout.CENTER)

        ret = vortex.showInDialog(mypanel, "ADME")
        if ret == vortex.OK:
            return True
        else:
            return False

    def parse_result(self, resultStr):
        return

    def calculateAll(self):
        this = self
        selectedProperties = []
        for key in self.checkboxDict.keys():
            if self.checkboxDict[key].isSelected():
                selectedProperties.append(key)
        if len(selectedProperties) == 0:
            return

        nrow = vtable.getRealRowCount()
        columnDict = {}
        for property in selectedProperties:
            column = vtable.findColumnWithName(property, True, this.columnTypeDict[property])
            columnDict[property] = column

        class ConvertTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setNote("Running ...")
                self.progressMonitor.setProgress(1)
                self.progressMonitor.setMillisToPopup(0)
                self.progressMonitor.setMillisToDecideToPopup(0)
                self.molUtilies = MolUtilities()

            def doInBackground(self):
                self.super__publish([(1,"Generating input ...")])
                mol_list = []
                for row in range(0, nrow):
                    if self.progressMonitor.isCanceled():
                        break
                    self.super__publish([(25*row/nrow,"process compound %d"%(row))])
                    try:
                        mfm = vtable.getMolFileManager()
                        mol_string = mfm.getMolFileAtRow(row)
                        mol = self.molUtilies.convertToMol(mol_string)
                        mol_list.append(mol)
                    except:
                        print >>sys.stderr, columnDict
                mol_input_string = self.molUtilies.convertMolsToSDF(mol_list)

                self.super__publish([(25,"contacting server ...")])
                #vortex.alert(" ".join(selectedProperties))
                # ["ALogD", "MoKa LogP", "MoKa LogD", "MoKa Basic pKa", "MoKa Acidic pKa",
                #  "ChemAxon LogP", "ChemAxon LogD", "ChemAxon Basic pKa",
                #  "Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)",
                #  "Cellarity_hERG_IC50(uM)"
                #  ]
                calculateALogD = "ALogD" in selectedProperties
                for mol_id, mol in enumerate(mol_list):
                    self.super__publish([(25 + 50 * mol_id / len(mol_list), "calculating property for mol %s" % (mol_id))])
                    propertyDict = {}
                    smilesDict = {}
                    smiles = self.molUtilies.molToSmiles(mol)
                    if smiles is not None:
                        if "\\" in smiles:
                            smiles = smiles.replace("\\","")
                        smilesDict[0] = smiles
                    if len(smilesDict) == 0:
                        continue

                    for idx,property in enumerate(selectedProperties):
                        if property == "ALogD":
                            continue
                        column = columnDict[property]
                        if property in ["ChemAxon LogP", "ChemAxon LogD", "ChemAxon Basic pKa", "ChemAxon Acidic pKa"]:
                            if property not in propertyDict:
                                acidic_pka,basic_pka = this.calculate_pKa_chemaxon(mol)
                                if acidic_pka is not None:
                                    propertyDict["ChemAxon Acidic pKa"] = acidic_pka[0]
                                else:
                                    propertyDict["ChemAxon Acidic pKa"] = None

                                if basic_pka is not None:
                                    propertyDict["ChemAxon Basic pKa"] = basic_pka[0]
                                else:
                                    propertyDict["ChemAxon Basic pKa"] = None

                                chemaxon_logd = this.calculate_cLogD74_chemaxon(mol)
                                propertyDict["ChemAxon LogD"] = chemaxon_logd
                                chemaxon_logp = this.calculate_cLogP_chemaxon(mol)
                                propertyDict["ChemAxon LogP"] = chemaxon_logp

                        if property in ["MoKa LogP", "MoKa LogD", "MoKa Basic pKa", "MoKa Acidic pKa"]:
                            if property not in propertyDict:
                                result_dict = this.calculate_MoKa(smilesDict)
                                #vortex.alert(str(result_dict))
                                if len(result_dict) > 0:
                                    basic_pka = result_dict['0']['most_basic']
                                    acidic_pka = result_dict['0']['most_acidic']
                                    logd74 = result_dict['0']['logd']
                                    log_sol = result_dict['0']['log_solubility']
                                    molar_sol = result_dict['0']['molar_solubility']
                                    propertyDict["MoKa LogP"] = 99.9
                                    propertyDict["MoKa LogD"] = logd74
                                    propertyDict["MoKa Basic pKa"] = basic_pka
                                    propertyDict["MoKa Acidic pKa"] = acidic_pka
                                    propertyDict["MoKa log Solubility"] = log_sol

                        if property == "Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)":
                            if property not in propertyDict:
                                result_dict = this.calculate_Perm(smilesDict)
                                if len(result_dict) > 0:
                                    perm = result_dict['0']['Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)']
                                    propertyDict['Cellarity Model PGP-KO PAPP AB (10^6 cm/sec)'] = perm

                        if property == "hepatocyte stability":
                            if property not in propertyDict:
                                result_dict = this.calculate_HepatocyteStability(smilesDict)
                                if len(result_dict) > 0:
                                   hep_stab = result_dict['0']['predicted_label']
                                   propertyDict[property] = hep_stab

                        if property == "Cellarity_hERG_IC50(uM)":
                            if property not in propertyDict:
                                result_dict = this.calculate_hERG(smilesDict)
                                if len(result_dict) > 0:
                                    hERG = result_dict['0']['Cellarity_hERG_IC50(uM)']
                                    propertyDict['Cellarity_hERG_IC50(uM)'] = hERG

                        if property in propertyDict:
                            try:
                                value = propertyDict[property]
                                if value is not None:
                                    if this.columnTypeDict[property] == vortex.DOUBLE:
                                        column.setValueFromString(mol_id, "%5.2f"%propertyDict[property])
                                    else:
                                        column.setValueFromString(mol_id,str(propertyDict[property]))
                            except:
                                pass
                if calculateALogD:
                    # f = open("/Users/jfeng/junk_junk.sdf","w")
                    # f.write(mol_input_string)
                    #print(mol_input_string,file=f)
                    column = columnDict["ALogD"]
                    result = this.xmlrpc_client.execute("alogd", [mol_input_string])
                    parser = JSONParser()
                    json_obj = parser.parse(result)
                    for row in range(0, nrow):
                        column.setValueFromString(row, "%f" % (json_obj.get("%d" % row)))

            def process(self, pairs):
                    pair = pairs[len(pairs)-1]
                    self.progressMonitor.setProgress(pair[0])
                    self.progressMonitor.setNote(pair[1])

            def done(self):
                self.progressMonitor.close()
                try:
                    self.get()  #raise exception if abnormal completion
                    vtable.fireTableStructureChanged()
                    vortex.alert("Finished.")
                except ExecutionException, e:
                    vortex.alert(e.getMessage())

        task = ConvertTask()
        task.execute()


if __name__ == "__main__":
    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        adme_runner = ADMERunner()
        if adme_runner.show_interface():
            adme_runner.calculateAll()