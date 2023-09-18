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
import java.io.ByteArrayInputStream as ByteArrayInputStream
import java.awt.event.ActionListener as ActionListener
import java.awt.event.KeyListener as KeyListener
import java.awt.event.KeyEvent as KeyEvent
import javax.swing.JOptionPane as JOptionPane

import java.beans.PropertyChangeListener as PropertyChangeListener

sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

import openeye.oechem as oechem
import openeye.oechem.OEGraphMol as OEGraphMol

import chemaxon.marvin.beans.MSketchPane as MSketchPane
import chemaxon.formats.MolImporter as MolImporter
import chemaxon.formats.MolExporter as MolExporter

import org.apache.xmlrpc.client.XmlRpcClient as XmlRpcClient
import org.apache.xmlrpc.client.XmlRpcClientConfigImpl as XmlRpcClientConfigImpl
import urllib, subprocess
import org.json.simple.JSONArray as JSONArray
import org.json.simple.JSONObject as JSONObject
import org.json.simple.parser.JSONParser as JSONParser
import org.json.simple.parser.ParseException as ParseException

class MolUtilities():
    def __init__(self):
        return

    def convertMolsToSDF(self,mols):
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat.SDF)
        ofs.openstring()
        for mol in mols:
            oemol = OEGraphMol(mol)
            oechem.oechem.OEWriteMolecule(ofs,oemol)
        return ofs.GetString()

    def convertToChemaxonMol(self,molString):
        mol = MolImporter.importMol(molString,"mol")
        return mol

    def convertToMol(self, molString):
        ifs = oechem.oemolistream()
        ifs.SetFormat(oechem.OEFormat.SDF)
        ifs.openstring(molString)
        mol = oechem.OEGraphMol()
        oechem.oechem.OEReadMolecule(ifs, mol)
        return mol

    def molToSmiles(self, mol):
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat.SMI)
        ofs.openstring()
        oechem.oechem.OEWriteMolecule(ofs, mol)
        return ofs.GetString()

    def molToSDF(self, mol):
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat.SDF)
        ofs.openstring()
        oechem.oechem.OEWriteMolecule(ofs, mol)
        return ofs.GetString()

    def get3DMol(self, mol):
        mol = oechem.OEMol(mol)
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat.SDF)
        ofs.openstring()
        if mol.GetDimension()!=3:
            omega = oeomega.OEOmega()
            omega.SetStrictStereo(False)
            omega.SetMaxConfs(1)
            if omega.call(mol):
                oechem.oechem.OEWriteMolecule(ofs, mol)
                return ofs.GetString()
            else:
                return None
        else:
            oechem.oechem.OEWriteMolecule(ofs, mol)
            return ofs.GetString()

    def getMolListFromSDFString(self, sdfString):
        oemols = []
        ifs = oechem.oemolistream()
        ifs.SetFormat(oechem.OEFormat.SDF)
        ifs.openstring(sdfString)
        mol = OEGraphMol()
        while oechem.oechem.OEReadMolecule(ifs, mol):
            oemols.append(OEGraphMol(mol))
        return oemols


class ADMERunner:
    def __init__(self):
        self.columnTypeDict = {
            "hydrogen-bond acceptors":Column.COLUMN_INT,
            "hydrogen-bond donors":Column.COLUMN_INT,
            "molecular weight":Column.COLUMN_DOUBLE,
            "2d PSA":Column.COLUMN_DOUBLE,
            "sum of formal charges":Column.COLUMN_DOUBLE,
            "number of rings":Column.COLUMN_INT,
            "No. Aromatic Rings":Column.COLUMN_INT,
            "No. Aromatic Ring Systems":Column.COLUMN_INT,
            "rotatable bonds":Column.COLUMN_INT,
            "Lipinski violations":Column.COLUMN_INT,
            "CLogP":Column.COLUMN_DOUBLE,
            "ChemAxon LogP":Column.COLUMN_DOUBLE,
            "ChemAxon LogD":Column.COLUMN_DOUBLE,
            "ChemAxon pKa":Column.COLUMN_DOUBLE,
            "hERG":Column.COLUMN_DOUBLE,
            "Papp A->B":Column.COLUMN_DOUBLE,
            "General Metabolic Stability: T1/2 (min) (Mouse)":Column.COLUMN_DOUBLE
        }
        self.properties = ["hydrogen-bond acceptors","hydrogen-bond donors","molecular weight","2d PSA",
                           "sum of formal charges","number of rings", "No. Aromatic Rings","No. Aromatic Ring Systems",
                           "rotatable bonds", "Lipinski violations",
                           "CLogP","ChemAxon LogP",
                           "ChemAxon LogD",
                           "ChemAxon pKa"]
#        self.adme_models = ["CNS MPO","hERG","Papp A->B","General Metabolic Stability: T1/2 (min) (Mouse)"]
        self.adme_models = ["hERG","Papp A->B","General Metabolic Stability: T1/2 (min) (Mouse)"]
        self.model_method_dict = {"hERG":"hERG","Papp A->B":"mdck","General Metabolic Stability: T1/2 (min) (Mouse)":"msm"}
        self.xmlrpc_client = XmlRpcClient()
        self.xmlrpc_config = XmlRpcClientConfigImpl()
        self.xmlrpc_config.setServerURL(URL("http://10.74.2.128:9555"))
        self.xmlrpc_client.setConfig(self.xmlrpc_config)
        self.checkboxDict = {}

    def build_property_panel(self,title,properties, ncolumn):
        nrow = len(properties)/2+1
        p = JPanel(GridLayout(nrow, ncolumn))
        p.setBorder(TitledBorder(title))
        for property in properties:
            cb = JCheckBox(property, True)
            p.add(cb)
            self.checkboxDict[property] = cb
        return p

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
        #p1 = self.build_property_panel("Properties",self.properties,3)
        p2 = self.build_property_panel("ADME",self.adme_models,2)
        mainPanel = JPanel()
        layout = BoxLayout(mainPanel,BoxLayout.Y_AXIS)
        mainPanel.setLayout(layout)
        #mainPanel.add(p1)
        mainPanel.add(p2)
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
            column = vtable.findColumnWithName(property, 1, this.columnTypeDict[property])
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
                        oechem.oechem.OEClearSDData(mol)
                        oechem.oechem.OESetSDData(mol,"moka_id","%d"%row)
                        mol_list.append(mol)
                    except:
                        print >>sys.stderr, columnDict
                mol_input_string = self.molUtilies.convertMolsToSDF(mol_list)

                self.super__publish([(25,"contacting server ...")])
                for idx,property in enumerate(selectedProperties):
                    self.super__publish([(25+50*idx/len(selectedProperties),"calculating property %s"%(property))])
                    column = columnDict[property]
                    result = this.xmlrpc_client.execute(this.model_method_dict[property],[mol_input_string])
                    parser = JSONParser()
                    json_obj = parser.parse(result)
                    for row in range(0,nrow):
                        column.setValueFromString(row,"%f"%(json_obj.get("%d"%row)))

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