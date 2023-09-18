__author__ = 'jfeng1'
import sys
import os
import traceback, tempfile, urllib2
#import json #only work under jython2.7.0
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
import com.dotmatics.vortex.VortexApplication as VortexApplication
import com.dotmatics.vortex.table.IntArray as IntArray
from com.dotmatics.chemlink import ChemLink
import com.dotmatics.vortex.table.IntArray as IntArray
import javax.swing.SwingUtilities as SwingUtilities
import java.lang.String as String
import org.python.core.util.StringUtil as StringUtil
import javax.swing.event.ChangeListener as ChangeListener
from javax.swing import BorderFactory, ProgressMonitor, SwingWorker, JScrollPane, JComboBox, JPanel, JTable, ImageIcon, \
    JLabel, JFileChooser, JCheckBox, JSlider
from javax.swing import JButton, JFrame, JTextField
from java.awt import GridBagLayout as awtGridBagLayout, BorderLayout, FlowLayout, Dimension
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
import java.util.Timer as Timer
import javax.swing.JDialog as JDialog
import javax.swing.JList as JList
import javax.swing.ListSelectionModel as ListSelectionModel
import javax.swing.event.ListSelectionListener as ListSelectionListener
import javax.swing.JScrollPane as JScrollPane
import javax.swing.JTextArea as JTextArea

sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

import org.json.simple.JSONObject as JSONObject
import org.json.simple.JSONArray as JSONArray
import org.json.simple.parser.JSONParser as JSONParser

import openeye.oechem as oechem
import openeye.oegraphsim as oegraphsim
import openeye.oedocking as oedocking
import openeye.oebio as oebio
import openeye.oeomega as oeomega
import org.apache.batik.transcoder.TranscoderInput as TranscoderInput
import org.apache.batik.transcoder.TranscodingHints as TranscodingHints
import org.apache.batik.transcoder.image.ImageTranscoder as ImageTranscoder
import openeye.oechem as oechem
import openeye.oedepict as oedepict
import openeye.oegrapheme as oegrapheme
import openeye.oeshape as oeshape
import openeye.oechem.OEGraphMol as OEGraphMol
import openeye.oedepict.OE2DMolDisplayOptions as OE2DMolDisplayOptions
import openeye.oeshape.OEColorForceField as OEColorForceField
import openeye.oeshape.OEColorFFType as OEColorFFType
import openeye.oegraphsim as oegraphsim
import openeye.oedepict.OE2DMolDisplay as OE2DMolDisplay
import openeye.oedepict.OEImageBase as OEImageBase
import openeye.oegrapheme.OE2DPropMap as OE2DPropMap
import openeye.oechem.OEColor as OEColor

import urllib, subprocess, xmlrpclib

PORT = 9527
MODELING_SERVER = "http://modeling.corp.biogen.com:%d/" % PORT

class PymolClientException(Exception):
    def __init__(self, value):
        self.value = value
        return

    def __str__(self):
        return repr(self.value)


class PymolClient():
    def __init__(self):
        self.status = 0
        self.url_prefix = "http://localhost:8080/apply/pymol.cmd"
        self.test_url = "http://localhost:8080/apply/pymol.cmd.pwd"
        try:
            urllib.urlopen(self.test_url)
        except:
            self.pymolExePath = vortex.getProperty("pymolexe")
            self.pymolDirectory = vortex.getProperty("pymolpath")
            if self.pymolExePath is None or self.pymolDirectory is None:
                self.pymolExePath, self.pymolDirectory = self.getInput()
            if self.pymolDirectory is not None and self.pymolDirectory is not None:
                vortex.setProperty("pymolexe", self.pymolExePath)
                vortex.setProperty("pymolpath", self.pymolDirectory)
                if Util.getPlatform() == Util.PlatformIsMac:
                    self.process = subprocess.Popen(["open", "-a", "%s" % self.pymolExePath,
                                                     "%s/pymol/modules/web/pymolhttpd.py" % self.pymolDirectory])
                elif Util.getPlatform() == Util.PlatformIsWindows:
                    self.process = subprocess.Popen(
                        ["%s" % self.pymolExePath, "%s\\modules\\web\\pymolhttpd.py" % self.pymolDirectory])
                else:
                    self.status = 1
                    vortex.alert("Platform not supported.")
                    return
                self.process.poll()
            else:
                self.status = 1
                vortex.alert("Failed to find pymol in the directory specified.")
        return

    def isSuccess(self):
        return self.status == 0

    def getInput(self):
        if Util.getPlatform() == Util.PlatformIsMac:
            vortex.alert("MacPymol not found, please browse to its location.")
        else:
            vortex.alert("pymol.exe not found, please browse to its location.")

        fileChooser = JFileChooser()
        fileChooser.setName("Select PyMol Location")
        val = fileChooser.showOpenDialog(None)
        if (val == 0):
            file = fileChooser.selectedFile
            if "pymol" not in file.getName().lower():
                return None, None
            if Util.getPlatform() == Util.PlatformIsMac:
                return file.getAbsolutePath(), file.getAbsolutePath()
            elif Util.getPlatform() == Util.PlatformIsWindows:
                return file.getAbsolutePath(), file.getParent()
            else:
                return None, None
        else:
            return None, None

    def open_url(self, url, time_out=1):
        success = False
        time_elapsed = 0
        while (success is False and time_elapsed <= time_out):
            try:
                urllib2.urlopen(url)
                success = True
            except:
                sleep(1)
                time_elapsed += 1
                success = False
        return success

    def open_file(self, filename, multiplex):
        filename1 = urllib.quote(filename)
        self.open_url("%s.reinitialize" % self.url_prefix)  # clear everything
        url = "%s.load?filename=%s&multiplex=%d" % (self.url_prefix, filename1, multiplex)
        self.open_url(url)
        self.open_url("%s.set?name=valence&value=1" % self.url_prefix)  # turn valence
        self.open_url("%s.hide?representation=lines&selection=all" % self.url_prefix)  # undisplay lines
        self.open_url("%s.show?representation=sticks" % self.url_prefix)  # just display sticks
        self.open_url("%s.hide?representation=sticks&selection=h." % self.url_prefix)  # undisplay hydrogen

    def enable(self, molname):
        molname1 = urllib.quote(molname)
        url = "%s.enable?name=%s" % (self.url_prefix, molname1)
        self.open_url(url)

    def disable(self, molname):
        molname1 = urllib.quote(molname)
        url = "%s.disable?name=%s" % (self.url_prefix, molname1)
        self.open_url(url)

    def show(self, representation, selectionName):
        selectionName1 = urllib.quote(selectionName)
        representation1 = urllib.quote(representation)
        url = "%s.show?representation=%s&selection=%s" % (self.url_prefix, representation1, selectionName1)
        self.open_url(url)

    def hide(self, representation, selectionName):
        selectionName1 = urllib.quote(selectionName)
        representation1 = urllib.quote(representation)
        url = "%s.hide?representation=%s&selection=%s" % (self.url_prefix, representation1, selectionName1)
        self.open_url(url)

    def zoom(self, selectionName):
        selectionName1 = urllib.quote(selectionName)
        url = "%s.zoom?selection=%s&state=-1" % (self.url_prefix, selectionName1)
        self.open_url(url)

def getAstexColor(rowId):
    colors = ["gray","yellow","magenta","cyan","pink","orange"]
    if rowId <0:
        return colors[0]
    else:
        return colors[rowId%len(colors)]


def getActiveSiteInteractionImage(receptor_mol, ligand_mol, width=1000, height=1000):
    active_site = oebio.OEFragmentNetwork(receptor_mol,ligand_mol)
    oedocking.oedocking.OEAddDockingInteractions(active_site)
    oegrapheme.oegrapheme.OEPrepareActiveSiteDepiction(active_site)
    image = oedepict.OEImage(width,height)
    center_frame = oedepict.OEImageFrame(image,width*0.8,height*0.8, oedepict.OE2DPoint(0.0,0.0))
    legend_frame = oedepict.OEImageFrame(image,width*0.2,width*0.2,oedepict.OE2DPoint(width*0.8,0.0))
    opts = oegrapheme.OE2DActiveSiteDisplayOptions(center_frame.GetWidth(),center_frame.GetHeight())
    active_site_display = oegrapheme.OE2DActiveSiteDisplay(active_site,opts)
    oegrapheme.oegrapheme.OERenderActiveSite(center_frame,active_site_display)
    legend_opts = oegrapheme.OE2DActiveSiteLegendDisplayOptions(10,1)
    oegrapheme.oegrapheme.OEDrawActiveSiteLegend(legend_frame,active_site_display,legend_opts)
    return image

class MyTranscoder(ImageTranscoder):
    def __init__(self):
        self.image = None

    def createImage(self, w, h):
        self.image = BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB)
        return self.image

    def getImage(self):
        return self.image

def getInput():
    fileChooser = JFileChooser()
    val = fileChooser.showSaveDialog(None)
    if val == JFileChooser.APPROVE_OPTION:
        f = fileChooser.selectedFile
        return str(f.getAbsolutePath())
    else:
        return None


import org.apache.batik.dom.svg.SVGDOMImplementation as SVGDOMImplementation
import org.apache.batik.util.SVGConstants as SVGConstants


class SVGTableCellRenderer(DefaultTableCellRenderer):
    def __init__(self):
        return

    def getTableCellRendererComponent(self, jtable, imageString, isSelected, hasFocus, row, column):
        if imageString is not None and len(imageString) > 0:
            transcoder = MyTranscoder()
            rect = jtable.getCellRect(row, column, True)
            svgImageArray = String(imageString).getBytes()
            transcoderInput = TranscoderInput(ByteArrayInputStream(svgImageArray))
            hints = TranscodingHints()
            hints.put(ImageTranscoder.KEY_WIDTH, java.lang.Float(rect.getWidth()))
            hints.put(ImageTranscoder.KEY_HEIGHT, java.lang.Float(rect.getHeight()))
            hints.put(ImageTranscoder.KEY_DOM_IMPLEMENTATION, SVGDOMImplementation.getDOMImplementation())
            hints.put(ImageTranscoder.KEY_DOCUMENT_ELEMENT_NAMESPACE_URI, SVGConstants.SVG_NAMESPACE_URI)
            hints.put(ImageTranscoder.KEY_DOCUMENT_ELEMENT, SVGConstants.SVG_SVG_TAG)
            hints.put(ImageTranscoder.KEY_XML_PARSER_VALIDATING, False)

            transcoder.setTranscodingHints(hints)
            transcoder.transcode(transcoderInput, None)
            imageIcon = ImageIcon(transcoder.getImage())
            self.setIcon(imageIcon)
        return self


class ModelingServer:
    def __init__(self):
        self.server =  xmlrpclib.Server(MODELING_SERVER)

    def get_docking_grids(self):
        grids = self.server.getAvailableDockingGrids()
        return grids

    def dock(self,grid_name, ligandMolString, numPoses):
        return self.server.docking(grid_name,ligandMolString,numPoses)

    def generatePymolSession(self,receptorStr,referenceStr, ligands, filename):
        dict = JSONObject()
        dict.put('receptor',receptorStr)
        dict.put('reference',referenceStr)
        ligList = JSONArray()
        for ligand in ligands:
            ligList.add(ligand)
        dict.put('ligands',ligList)
        res = self.server.generatePymolSession(dict.toJSONString())
        if res is not None:
            f = open(filename, "wb")
            f.write(res.data)
            f.close()
            return True
        return False


class MolUtilities():
    def __init__(self):
        return

    def convertToMol(self, mol_string, mol_format=oechem.OEFormat.SDF):
        ifs = oechem.oemolistream()
        ifs.SetFormat(mol_format)
        ifs.openstring(mol_string)
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

    def get2DMol(self,mol):
        mol = oechem.OEMol(mol)
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat.SDF)
        ofs.openstring()
        if mol.GetDimension()!=2:
            if oedepict.oedepict.OEPrepareDepiction(mol,True,True): #clear coordinates and suppress hydrogen
                oechem.oechem.OEWriteMolecule(ofs, mol)
                return ofs.GetString()
            else:
                return None
        else:
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


class JobRunner():
    def __init__(self):
        self.progressMonitor = ProgressMonitor(vortex, "Progress", "Progress", 0, 100)
        self.progressMonitor.setMillisToPopup(0)
        self.progressMonitor.setMillisToDecideToPopup(0)
        self.isCancelled = False
        self.selectedColumn = None
        self.numPoses = 100
        self.nameCB = None
        self.modelingServer = ModelingServer()
        self.inputMol = None
        self.timer = Timer()
        return

    def show_interface(self):

        dialog = JDialog(VortexApplication.vortexApplication,"Select Receptor For Docking",True)
        dialog.setSize(Dimension(800,600))
        dialog.setLocationRelativeTo(None)
        dialog.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE)

        btnPanel = JPanel()
        okButton = JButton("Ok")
        this = self
        class OkListener(ActionListener):
            def actionPerformed(self,e):
                this.isCancelled = False
                dialog.setVisible(False)
        okButton.addActionListener(OkListener())

        cancelButton = JButton("Cancel")
        class CancelListener(ActionListener):
            def actionPerformed(self,e):
                this.isCancelled = True
                dialog.setVisible(False)
        cancelButton.addActionListener(CancelListener())
        okButton.setEnabled(False)

        btnPanel.add(okButton)
        btnPanel.add(cancelButton)

        grids = self.modelingServer.get_docking_grids()
        gridList = JList(grids)
        gridList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION)
        class MyListSelectionListener(ListSelectionListener):
            def valueChanged(self,e):
                if not e.getValueIsAdjusting():
                    if gridList.getSelectedIndex == -1:
                        okButton.setEnabled(False)
                    else:
                        okButton.setEnabled(True)
        gridList.getSelectionModel().addListSelectionListener(MyListSelectionListener())
        mainPanel = JPanel(BorderLayout())
        mainPanel.add(JScrollPane(gridList),BorderLayout.CENTER)
        mainPanel.add(btnPanel,BorderLayout.SOUTH)
        dialog.getContentPane().add(mainPanel)
        dialog.setVisible(True)
        if self.isCancelled:
            return False

        idx = gridList.getSelectedIndex()
        if idx == -1:
            return False

        self.grid_name = grids[idx]

        chemLink = ChemLink()
        if vtable is not None:
            selectedRows = vtable.getSelectedRows()
            if selectedRows.size() >= 1:
                selectedRow = selectedRows.get(0)
                mfm = vtable.getMolFileManager()
                self.inputMol = mfm.getMolFileAtRow(selectedRow)
                if self.inputMol is not None:
                    mol_utilities = MolUtilities()
                    mol = mol_utilities.convertToMol(self.inputMol)
                    molString = mol_utilities.get2DMol(mol)
                    chemLink.setMolFile(molString)
                    class MyRunnable(java.lang.Runnable):
                        def run(self):
                            chemLink.setMolFile(molString)
                            chemLink.repaint()
                    SwingUtilities.invokeLater(MyRunnable())

        chemLink.setPreferredSize(Dimension(600, 400))
        p = JPanel(BorderLayout())
        p.add(chemLink, BorderLayout.CENTER)
        p2 = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        self.slider = JSlider(0,10)
        self.slider.setMajorTickSpacing(1)
        self.slider.setPaintTicks(True)
        self.slider.setPaintLabels(True)
        self.slider.setLabelTable(self.slider.createStandardLabels(5))
        textField = JTextField(5)
        textField.setText("5")
        textField.setEditable(False)
        slider1 = self.slider
        class SliderChangeListener(ChangeListener):
            def stateChanged(self, e):
                textField.setText("%d"%slider1.getValue())

        self.slider.addChangeListener(SliderChangeListener())
        p2.add(JLabel("No. of Poses:"))
        p2.add(self.slider)
        p2.add(textField)
        p.add(p2,BorderLayout.SOUTH)
        p.setBorder(BorderFactory.createTitledBorder("Click to edit structure:"))

        ret = vortex.showInDialog(p, "Sketch input molecule for Docking")
        if ret == 0:
            if chemLink.getMolFile() is None or len(chemLink.getMolFile().strip()) == 0:
                return False
            else:
                self.inputMol = chemLink.getMolFile()
                return True
        else:
            return False

    def run(self):
        this = self
        this.numPoses = this.slider.getValue()
        if this.numPoses == 0:
            this.numPoses = 5

        class MyTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self, progressMonitor):
                self.progressMonitor = progressMonitor
                self.molFileManager = {}
                self.molUtilities = MolUtilities()
                self.myProgress = 1
                self.receptorPdb = None
                self.referenceSdf = None
                SwingWorker.__init__(self)

            def doInBackground(self):
                self.super__publish([("Docking molecule ...", self.myProgress)])
                try:
                    rs = this.modelingServer.dock(this.grid_name,this.inputMol, this.numPoses)
                    parser = JSONParser()
                    result = parser.parse(rs)
                    print >> sys.stderr, rs
                    ligandSdf = result['ligand']
                    self.receptorPdb = result['receptor']
                    receptorMol = self.molUtilities.convertToMol(self.receptorPdb,oechem.OEFormat.PDB)
                    self.referenceSdf = result['reference']

                    vortexTableModel = VortexTableModel()
                    structureColumn = vortexTableModel.findColumnWithName("MOLFILE", 1, Column.COLUMN_STRING)
                    imageColumn = vortexTableModel.findColumnWithName("Image", 1, Column.COLUMN_STRING)
                    nameColumn = vortexTableModel.findColumnWithName("Name", 1, Column.COLUMN_STRING)
                    scoreColumn = vortexTableModel.findColumnWithName("Glide Score", 1,Column.COLUMN_DOUBLE)
                    refMol3d = self.molUtilities.convertToMol(self.referenceSdf)
                    self.myProgress = self.myProgress + 10
                    self.super__publish([("Generating image ...", self.myProgress)])
                    refMol2d = OEGraphMol(refMol3d)
                    oegrapheme.oegrapheme.OEPrepareDepictionFrom3D(refMol2d)
                    images = getActiveSiteInteractionImage(receptorMol,refMol3d)
                    # mol2d, images = getOverlayImagesRef(refMol3d)
                    self.myProgress = self.myProgress + 30
                    self.super__publish([("Finished generating image ...", self.myProgress)])

                    structureColumn.setValueFromString(0, self.molUtilities.molToSDF(refMol2d))
                    nameColumn.setValueFromString(0, refMol3d.GetTitle())
                    scoreColumn.setValueFromString(0,None)
                    array = oedepict.oedepict.OEWriteImageToByteArray("svg", images)
                    imageColumn.setValueFromString(0, String(array))
                    self.molFileManager[0] = self.referenceSdf
                    # fingerprint_ref = oegraphsim.OEFingerPrint()
                    # oegraphsim.oegraphsim.OEMakeMACCS166FP(fingerprint_ref, refMol3d)
                    # #                molFileManager.addMoleFile(String(this.refSDF))
                    drugs = self.molUtilities.getMolListFromSDFString(ligandSdf)
                    print >> sys.stderr,"Number of poses:%d"%len(drugs)
                    # for id, oe_mol in enumerate(drugs,start=1): # only works for jython 2.7
                    for id, oe_mol in enumerate(drugs):
                        id2 = id+1
                        self.super__publish([("building table ...", 10 + 90 * id / len(drugs))])
                        # fingerprint_mol = oegraphsim.OEFingerPrint()
                        # oegraphsim.oegraphsim.OEMakeMACCS166FP(fingerprint_mol, oe_mol)
                        # similarity = oegraphsim.oegraphsim.OETanimoto(fingerprint_ref, fingerprint_mol)

                        self.molFileManager[id2] = self.molUtilities.molToSDF(oe_mol)
                        mol2d = OEGraphMol(oe_mol)
                        oegrapheme.oegrapheme.OEPrepareAlignedDepictionFrom3D(mol2d, oe_mol, refMol2d, refMol3d)
                        overlay_images = getActiveSiteInteractionImage(receptorMol,oe_mol)
                        structureColumn.setValueFromString(id2, self.molUtilities.molToSDF(mol2d))
                        nameColumn.setValueFromString(id2, oe_mol.GetTitle())
                        array = oedepict.oedepict.OEWriteImageToByteArray("svg", overlay_images)
                        imageColumn.setValueFromString(id2, String(array))
                        scoreColumn.setValueFromString(id2,"%f"%round(float(oechem.oechem.OEGetSDData(oe_mol,"r_i_glide_gscore")),2))

                    return vortexTableModel

                except:
                    print "My Exception:"
                    traceback.print_exc()

                return None

            def process(self, pairs):
                pair = pairs[len(pairs) - 1]
                self.progressMonitor.setNote(pair[0])
                self.progressMonitor.setProgress(pair[1])

            def done(self):
                self.progressMonitor.close()
                try:
                    vortexTableModel = self.get()
                    if vortexTableModel is None:
                        raise ExecutionException(java.lang.Throwable("No result found."))
                    mainPanel = JPanel(BorderLayout())

                    import astex.MoleculeViewer as MoleculeViewer
                    import astex.FILE as AstexFile
                    import astex.MoleculeIO as MoleculeIO

                    molViewer = MoleculeViewer()
                    molViewer.execute("set symmetry off;")
                    viewPanel = JPanel(BorderLayout())
                    viewPanel.add(molViewer, BorderLayout.CENTER)

                    scriptArea = JTextArea()
                    scriptArea.setLineWrap(True)
                    scriptArea.setWrapStyleWord(True)
                    scriptBtn = JButton("Script...")
                    dialog = JDialog(VortexApplication.vortexApplication,"Execute ",True)
                    dialog.setSize(Dimension(300,200))
                    dialog.setLocationRelativeTo(None)
                    dialog.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE)
                    dialog_btnPanel = JPanel()
                    okButton = JButton("Ok")
                    class OkListener(ActionListener):
                        def actionPerformed(self,e):
                            try:
                                molViewer.execute(scriptArea.getText())
                                molViewer.getMoleculeRenderer().repaint()
                                dialog.setVisible(False)
                            except:
                                vortex.alert("Wrong script!")
                    okButton.addActionListener(OkListener())

                    cancelButton = JButton("Cancel")
                    class CancelListener(ActionListener):
                        def actionPerformed(self,e):
                            dialog.setVisible(False)
                    cancelButton.addActionListener(CancelListener())
                    okButton.setEnabled(False)

                    dialog_btnPanel.add(okButton)
                    dialog_btnPanel.add(cancelButton)

                    dialog_main_panel = JPanel(BorderLayout())
                    dialog_main_panel.add(JScrollPane(scriptArea),BorderLayout.CENTER)
                    dialog_main_panel.add(dialog_btnPanel,BorderLayout.SOUTH)
                    dialog.getContentPane().add(dialog_main_panel)

                    class ScriptBtnActionListener(ActionListener):
                        def actionPerformed(self, e):
                            dialog.setVisible(True)

                    scriptBtn.addActionListener(ScriptBtnActionListener())

                    class ViewerKeyListener(KeyListener):
                        def keyTyped(self, e):
                            dialog.getRootPane().setDefaultButton(okButton)
                            okButton.setEnabled(True)

                        def keyPressed(self, e):
                            if e.getKeyCode() == KeyEvent.VK_ENTER:
                                okButton.doClick()
                            return

                        def keyReleased(self, e):
                            return

                    scriptArea.addKeyListener(ViewerKeyListener())

                    hbondCB = JCheckBox("H Bond", True)
                    class HBondActionListener(ActionListener):
                        def actionPerformed(self,e):
                            if hbondCB.isSelected():
                                molViewer.execute("distance -delete hbonds; define lig not (aminoacid or solvent); distance -name hbonds -from { group lig and (element 7 8)  or (solvent and sphere 4 around group lig)} -to {element 7 8 and (aminoacid or solvent)} -contact 0.3 -mode nbpairs -colour yellow -on 0.3 -off 0.4 -radius 0.04;")
                            else:
                                molViewer.execute("distance -delete hbonds;")
                            molViewer.getMoleculeRenderer().repaint()
                    hbondCB.addActionListener(HBondActionListener())

                    autoZoomCB = JCheckBox("Auto Zoom", True)

                    viewerBtnPanel = JPanel(FlowLayout(FlowLayout.CENTER, 20, 10))
                    viewerBtnPanel.add(hbondCB)
                    viewerBtnPanel.add(autoZoomCB)
                    viewerBtnPanel.add(scriptBtn)

                    viewPanel.add(viewerBtnPanel,BorderLayout.PAGE_END) #for testing script of astex viewer

                    visibleMolDict = {}
                    molDict = {}
                    molFileManager = self.molFileManager
                    molUtilities = self.molUtilities
                    pdbString = AstexFile(ByteArrayInputStream(String(self.receptorPdb).getBytes()))
                    receptor = MoleculeIO.readPDB(pdbString)
                    receptor.setName("receptor")
                    receptor.setDisplayed(1)
                    molViewer.addMolecule(receptor)
                    molViewer.execute("display lines off molecule receptor;")
                    molViewer.execute("display lines on molecule receptor and not element 1;")
                    molViewer.getMoleculeRenderer().repaint()
                    class FastRocsTable(VortexJTable2):
                        def vortexTableChanged(self, tableModel, action, i1, i2, o):
                            if action == "SelectionChanged" and viewPanel.isShowing():
                                rowids = vortexTableModel.getSelectedRows().getArray()
                                for rowId in rowids:
                                    print >> sys.stderr, "%s is selected."%rowId
                                    if not visibleMolDict.has_key(rowId):
                                        if not molDict.has_key(rowId):
                                            molFile = String(molFileManager[rowId])
                                            oemol = molUtilities.convertToMol(molFile)
                                            oechem.oechem.OESuppressHydrogens(oemol, True)
                                            molString = AstexFile(ByteArrayInputStream(String(molUtilities.molToSDF(oemol)).getBytes()))
                                            molecule = MoleculeIO.readMDLMol(molString)
                                            molDict[rowId] = molecule
                                            molecule.setName("molecule_at_row_%d" % rowId)
                                            print >> sys.stderr, "DEBUG:"+molecule.getName()
                                            molViewer.addMolecule(molecule)
                                        else:
                                            molecule = molDict[rowId]

                                        molecule.setDisplayed(1)
                                        visibleMolDict[rowId] = True
                                        molViewer.execute("display lines off molecule %s;"%molecule.getName())
                                        molViewer.execute("display cylinders on molecule %s;"%molecule.getName())
                                        molViewer.execute("color %s molecule %s and element 6;"%(getAstexColor(rowId),molecule.getName()))
                                        if autoZoomCB.isSelected():
                                            molViewer.execute("center molecule %s;"%molecule.getName())
                                        #                                        molRenderer = molViewer.getMoleculeRenderer()
                                        #                                        atomIterator = molRenderer.getAtomIterator()
                                        #                                        while atomIterator.hasMoreElements():
                                        #                                            atom = atomIterator.getNextAtom()
                                        #                                            if atom.getElement()==1 and atom.getBondCount()==1:
                                        #                                                nbr = atom.getBondedAtom(0)
                                        #                                                if nbr.getElement() == 6:
                                        #                                                    atom.setSelected(True)
                                        #                                        molViewer.execute("display cylinders off current;")

                                for rowId in visibleMolDict.keys():
                                    if not vortexTableModel.getRowSelected(rowId):
                                        del (visibleMolDict[rowId])
                                        if molDict.has_key(rowId):
                                            mr = molViewer.getMoleculeRenderer()
                                            mol = mr.getMolecule(molDict[rowId].getName())
                                            if mol is not None:
                                                mol.setDisplayed(0)
                                            # molViewer.execute("molecule remove %s;" % molDict[rowId].getName())
                                            # del (molDict[rowId])

                                if hbondCB.isSelected():
                                    molViewer.execute("distance -delete hbonds; define lig not (aminoacid or solvent); distance -name hbonds -from { group lig and (element 7 8)  or (solvent and sphere 4 around group lig)} -to {element 7 8 and (aminoacid or solvent)} -contact 0.3 -mode nbpairs -colour yellow -on 0.3 -off 0.4 -radius 0.04;")
                                else:
                                    molViewer.execute("distance -delete hbonds;")

                                molViewer.getMoleculeRenderer().repaint()

                    table = FastRocsTable()
                    table.setModel(vortexTableModel)

                    svgCellIdx = 2
                    column = table.getColumnModel().getColumn(svgCellIdx)
                    column.setCellRenderer(SVGTableCellRenderer())
                    mainPanel.add(JScrollPane(table), BorderLayout.CENTER)
                    btnPanel = JPanel(FlowLayout(FlowLayout.CENTER, 20, 10))
                    btn = JButton("View Selected In 3D")
                    btnPanel.add(btn)

                    btn4 = JButton("Export to SDF")
                    btnPanel.add(btn4)

                    btn5 = JButton("Export to pymol session")
                    btnPanel.add(btn5)

                    mainPanel.add(btnPanel, BorderLayout.PAGE_END)
                    ws = vortex.createWorkspace(vortexTableModel,
                                                "%s_%s" % ("Docking Result", vortex.getWorkspaceCount()))
                    ws.addComponent(mainPanel, "Docking")
                    ws.addComponent(viewPanel, "3D Viewer")
                    ws.updateLayout(VortexWorkspace.LAYOUT_GRID)
                    if vortexTableModel.getRowCount()>0:
                        vortexTableModel.setRowSelected(0,True)
                    if vortexTableModel.getRowCount()>1:
                        vortexTableModel.setRowSelected(1,True)
                    vortexTableModel.fireTableSelectionChanged()

                    class View3DActionListener(ActionListener):
                        def actionPerformed(self, e):
                            if not viewPanel.isShowing():
                                ws.addComponent(viewPanel, "3D Viewer")
                                ws.updateLayout(VortexWorkspace.LAYOUT_GRID)
                                table.vortexTableChanged(vortexTableModel, "SelectionChanged", 0, 0, None)
                                #                            vortex.showInDialog(viewPanel,"molecules")

                    btn.addActionListener(View3DActionListener())

                    class ExportSDFListener(ActionListener):
                        def actionPerformed(self, e):
                            rowids = vortexTableModel.getSelectedRows().getArray()
                            if len(rowids) == 0:
                                vortex.alert("No molecule selected.")
                            else:
                                filename = getInput()
                                if filename is not None:
                                    f = open(filename, "w")
                                    for rowId in rowids:
                                        molFile = molFileManager[rowId]
                                        f.write(molFile)
                                    f.close()
                                    vortex.alert("File written to %s" % filename)
                                else:
                                    vortex.alert("No file specified.")
                    btn4.addActionListener(ExportSDFListener())

                    receptorStr = self.receptorPdb
                    referenceStr = self.referenceSdf

                    class ExportPymolListener(ActionListener):
                        def actionPerformed(self,e):
                            rowids = vortexTableModel.getSelectedRows().getArray()
                            if len(rowids) == 0:
                                vortex.alert("No molecule selected.")
                            else:
                                filename = getInput()
                                if filename is not None:
                                    if not filename.endswith(".pse"):
                                        filename = "%s.pse"%filename

                                    ligandList = []
                                    for rowId in rowids:
                                        molFile = molFileManager[rowId]
                                        ligandList.append(molFile)
                                    if this.modelingServer.generatePymolSession(receptorStr,referenceStr,ligandList,filename):
                                        vortex.alert("File written to %s" % filename)
                                    else:
                                        vortex.alert("Failed to generate pymol session.")
                                else:
                                    vortex.alert("No file specified.")
                    btn5.addActionListener(ExportPymolListener())


                except ExecutionException, e:
                    traceback.print_exc()
                finally:
                    self.progressMonitor.close()
                return

        task = MyTask(self.progressMonitor)
        task.execute()

import javax.swing.JPopupMenu as JPopupMenu

JPopupMenu.setDefaultLightWeightPopupEnabled(False)
jobRunner = JobRunner()
if jobRunner.show_interface():
    jobRunner.run()
