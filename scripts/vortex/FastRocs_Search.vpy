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
import javax.swing.JOptionPane as JOptionPane

import java.beans.PropertyChangeListener as PropertyChangeListener

sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

import openeye.oechem as oechem
import openeye.oegraphsim as oegraphsim
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

import chemaxon.marvin.beans.MSketchPane as MSketchPane
import chemaxon.formats.MolImporter as MolImporter
import chemaxon.formats.MolExporter as MolExporter

import org.apache.xmlrpc.client.XmlRpcClient as XmlRpcClient
import org.apache.xmlrpc.client.XmlRpcClientConfigImpl as XmlRpcClientConfigImpl
import urllib, subprocess

#FASTROCS_SERVER = "http://rd3lmxp07:%d/RPC2" % PORT
PORT1 = 9527
PORT2 = 9528
PORT3 = 9529
FASTROCS_SERVER_ALL = "http://fragollino.biogen.com:%d/RPC2" % PORT1
FASTROCS_SERVER_EVOTEC = "http://fragollino.biogen.com:%d/RPC2" % PORT2
FASTROCS_SERVER_SOFTFOCUSDIVERSE = "http://fragollino.biogen.com:%d/RPC2" % PORT3


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


class OEOverlayView:
    def __init__(self, opts):
        self.refColorMol2D = OEGraphMol()
        self.fitColorHits = OEGraphMol()
        self.fitColorMisses = OEGraphMol()

        self.options = opts
        self.maxSelfOverlap = 0.0

    def setMols(self, ref3D, fit3D):
        if ref3D.GetDimension() != 3:
            oechem.OEThrow.Warning("OEOverlayView ref mol must have 3D coords")
            return False

        if fit3D.GetDimension() != 3:
            oechem.OEThrow.Warning("OEOverlayView fit mol must have 3D coords")
            return False
        self.refmol = OEGraphMol(ref3D)
        self.fitmol = OEGraphMol(fit3D)

        cff = OEColorForceField()
        cff.Init(OEColorFFType.ImplicitMillsDean)
        oeshape.oeshape.OEAddColorAtoms(self.refmol, cff)
        oeshape.oeshape.OECompressColorAtoms(self.refmol)

        # set up refmol2D
        self.refmol2D = OEGraphMol(self.refmol)
        oegrapheme.oegrapheme.OEPrepareDepictionFrom3D(self.refmol2D, True)
        # self.AddAtomArcs(self.refmol2D)

        # set up fitmol2D
        self.fitmol2D = OEGraphMol(self.fitmol)
        # generate fit coords aligned to ref based on 3D alignment
        print "Before Conversion"
        print MolUtilities().molToSDF(self.fitmol2D)
        oegrapheme.oegrapheme.OEPrepareAlignedDepictionFrom3D(self.fitmol2D, self.fitmol, self.refmol2D, self.refmol)
        print "After Conversion"
        print MolUtilities().molToSDF(self.fitmol2D)

        oeshape.oeshape.OEAddColorAtoms(self.fitmol, cff)
        oeshape.oeshape.OECompressColorAtoms(self.fitmol)

        oeshape.oeshape.OEAddColorAtoms(self.refmol2D, cff)
        oeshape.oeshape.OECompressColorAtoms(self.refmol2D)

        # set up color hit and miss molecules
        self.SetupColorHitsAndMisses()
        return True

    def UpdateOptions(self, img):
        if not self.refmol2D.IsValid() or not self.fitmol2D.IsValid():
            return

        self.options.SetWidth(img.GetWidth())
        self.options.SetHeight(img.GetHeight())
        self.options.SetScale(oedepict.OEScale.AutoScale)

        refScale = oegrapheme.oegrapheme.OEGetMoleculeSurfaceScale(self.refmol2D, self.options)
        fitScale = oedepict.oedepict.OEGetMoleculeScale(self.fitmol2D, self.options)
        self.options.SetScale(Math.min(fitScale, refScale))

    def DrawRefDepiction(self, img):
        if not self.refmol2D.IsValid():
            return

        self.UpdateOptions(img)
        disp = oedepict.OE2DMolDisplay(self.refmol2D, self.options)
        oedepict.oedepict.OERenderMolecule(img, disp)

    def DrawRefColorDepiction(self, img):
        if not self.refmol2D.IsValid():
            return

        self.UpdateOptions(img)
        disp = oedepict.OE2DMolDisplay(self.refmol2D, self.options)
        refColorDisp = oedepict.OE2DMolDisplay(self.refColorMol2D, self.options)
        oedepict.oedepict.OEOffsetMolDisplay(refColorDisp, disp)

        oegrapheme.oegrapheme.OEDraw2DSurface(disp)
        self.DrawColorAtoms(disp.GetLayer(oedepict.OELayerPosition.Below), refColorDisp, oedepict.OEFill.On)
        oedepict.oedepict.OERenderMolecule(img, disp)

    def DrawFitShapeOverlay(self, img):
        if not self.refmol2D.IsValid() or not self.fitmol2D.IsValid():
            return

        self.UpdateOptions(img)
        refDisp = oedepict.OE2DMolDisplay(self.refmol2D, self.options)
        fitDisp = oedepict.OE2DMolDisplay(self.fitmol2D, self.options)
        oedepict.oedepict.OEOffsetMolDisplay(fitDisp, refDisp)

        propMap = oegrapheme.OE2DPropMap()
        propMap.SetMinValue(-1.0)
        propMap.SetMaxValue(self.maxSelfOverlap)
        propMap.SetPositiveColor(oechem.oechem.getOEBrown())
        propMap.SetLegendLocation(oegrapheme.OELegendLocation.Hidden)
        propMap.Render(fitDisp, "Overlap")

        oegrapheme.oegrapheme.OEDraw2DSurface(fitDisp, refDisp)
        oedepict.oedepict.OERenderMolecule(img, fitDisp)

    def DrawFitColorOverlay(self, img):
        if not self.refmol2D.IsValid() or not self.fitmol2D.IsValid():
            return

        self.UpdateOptions(img)
        refDisp = oedepict.OE2DMolDisplay(self.refmol2D, self.options)
        fitDisp = oedepict.OE2DMolDisplay(self.fitmol2D, self.options)
        oedepict.oedepict.OEOffsetMolDisplay(fitDisp, refDisp)

        oegrapheme.oegrapheme.OEDraw2DSurface(fitDisp, refDisp)

        # hits
        print >> sys.stderr, "hit:%d" % self.fitColorHits.NumAtoms()
        if self.fitColorHits.NumAtoms() > 0:
            hitDisp = oedepict.OE2DMolDisplay(self.fitColorHits, self.options)
            oedepict.oedepict.OEOffsetMolDisplay(hitDisp, refDisp)
            self.DrawColorAtoms(fitDisp.GetLayer(oedepict.OELayerPosition.Below), hitDisp, oedepict.OEFill.On)

        # misses
        print >> sys.stderr, "miss:%d" % self.fitColorMisses.NumAtoms()
        if self.fitColorMisses.NumAtoms() > 0:
            missDisp = oedepict.OE2DMolDisplay(self.fitColorMisses, self.options)
            oedepict.oedepict.OEOffsetMolDisplay(missDisp, refDisp)
            self.DrawColorAtoms(fitDisp.GetLayer(oedepict.OELayerPosition.Below), missDisp, oedepict.OEFill.Off)

        oedepict.oedepict.OERenderMolecule(img, fitDisp)

    def GetColorAtomColor(self, ctype):
        c = oechem.OEColor(oechem.oechem.getOEBlack())
        if ctype == 0:
            c = oechem.OEColor(oechem.oechem.getOEBrown())
        elif ctype == 1:
            c = oechem.OEColor(128, 128, 255)
        elif ctype == 2:
            c = oechem.OEColor(255, 89, 89)
        elif ctype == 3:
            c = oechem.OEColor(0, 0, 255)
        elif ctype == 4:
            c = oechem.OEColor(255, 0, 0)
        elif ctype == 5:
            c = oechem.OEColor(0, 255, 0)
        elif ctype == 6:
            c = oechem.OEColor(255, 204, 51)

        c.SetA(128)
        return c

    def DrawColorAtoms(self, image, disp, fill):
        colorRad = disp.GetOptions().GetAtomLabelFont().GetSize()
        pen = oedepict.OEPen()
        pen.SetFill(fill)
        for adisp in disp.GetAtomDisplays():
            ctr = adisp.GetCoords()
            colorType = oeshape.oeshape.OEGetColorType(adisp.GetAtom())
            c = self.GetColorAtomColor(colorType)
            pen.SetForeColor(c)
            pen.SetBackColor(c)
            image.DrawCircle(ctr, colorRad, pen)

    def AddAtomArcs(self, mol):
        for atom in mol.GetAtoms():
            b = oechem.oechem.getOEOliveGrey()
            p = oedepict.OEPen(b, b)
            p.SetLineStipple(oedepict.OEStipple.Dot)
            p.SetLineWidth(1.0)
            oegrapheme.oegrapheme.OESetSurfaceArcFxn(mol, atom, oegrapheme.OEDefaultArcFxn(p))

    def debug(self):
        refAtom3DIter = self.refmol.GetAtoms(oeshape.OEIsColorAtomPred())
        refAtom2DIter = self.refmol2D.GetAtoms(oeshape.OEIsColorAtomPred())
        print refAtom3DIter.IsValid(), refAtom2DIter.IsValid()

    def SetupColorHitsAndMisses(self):
        #       make a molecule of just the ref color atoms
        #        cff = OEColorForceField()
        #        cff.Init(OEColorFFType.ImplicitMillsDean)
        #        oeshape.oeshape.OEAddColorAtoms(self.refmol2D, cff)
        #        oeshape.oeshape.OECompressColorAtoms(self.refmol2D)
        #
        oeshape.oeshape.OEUncompressColorAtoms(self.refmol2D)
        self.refColorMol2D.Clear()
        for catom in self.refmol2D.GetAtoms(oeshape.OEIsColorAtomPred()):
            self.refColorMol2D.NewAtom(catom)
        self.refColorMol2D.SetDimension(2)
        oeshape.oeshape.OEUncompressColorAtoms(self.refmol)
        oeshape.oeshape.OEUncompressColorAtoms(self.fitmol)

        #        // make a molecule of hits and misses
        self.fitColorHits.Clear()
        self.fitColorHits.SetDimension(2)
        self.fitColorMisses.Clear()
        self.fitColorMisses.SetDimension(2)

        refAtom3DIter = self.refmol.GetAtoms(oeshape.OEIsColorAtomPred())
        refAtom2DIter = self.refmol2D.GetAtoms(oeshape.OEIsColorAtomPred())

        while refAtom3DIter.IsValid() and refAtom2DIter.IsValid():
            refAtom3D = refAtom3DIter.Target()
            refAtom2D = refAtom2DIter.Target()
            match = False
            for fitAtom3D in self.fitmol.GetAtoms(oeshape.OEIsColorAtomPred()):
                if oeshape.oeshape.OEGetColorType(refAtom3D) == oeshape.oeshape.OEGetColorType(fitAtom3D):
                    dist = oechem.oechem.OEGetDistance2(self.refmol, refAtom3D, self.fitmol, fitAtom3D)
                    if dist < 2.0:
                        self.fitColorHits.NewAtom(refAtom2D)
                        match = True
                        break

            if not match:
                self.fitColorMisses.NewAtom(refAtom2D)

            refAtom3DIter.Increment()
            refAtom2DIter.Increment()

        oeshape.oeshape.OECompressColorAtoms(self.refmol)
        oeshape.oeshape.OECompressColorAtoms(self.refmol2D)
        oeshape.oeshape.OECompressColorAtoms(self.fitmol)

    def AttachShapeOverlapData(self):
        maxSelfOverlap = 0.0

        #        // get shape match into propmap
        ov = oeshape.OEOverlap();
        ov.SetRefMol(self.refmol)

        #        // calculate max self overlap for an atom for scale
        ovres = oeshape.OEOverlapResults()
        refArray = oechem.OEFloatArray(self.refmol.GetMaxAtomIdx())
        ov.Overlap(self.refmol, ovres, refArray)
        for atom in self.refmol.GetAtoms():
            if refArray.getItem(atom.GetIdx()) > maxSelfOverlap:
                maxSelfOverlap = refArray.getItem(atom.GetIdx())

                #        // now calculate the atom overlap for the fit mol
        fitArray = oechem.OEFloatArray(self.fitmol.GetMaxAtomIdx())
        ov.Overlap(self.fitmol, ovres, fitArray)
        for atom in self.fitmol.GetAtoms():
            atom.SetDoubleData("Overlap", fitArray.getItem(atom.GetIdx()))


def getOverlayImages(refMol, fitMol):
    opts = OE2DMolDisplayOptions()
    overlayView = OEOverlayView(opts)
    overlayView.setMols(refMol, fitMol)
    img = oedepict.OEImage(300, 300)
    overlayView.DrawFitColorOverlay(img)
    return overlayView.fitmol2D, img


def getOverlayImagesRef(refMol):
    opts = OE2DMolDisplayOptions()
    overlayView = OEOverlayView(opts)
    overlayView.setMols(refMol, refMol)
    img = oedepict.OEImage(300, 300)
    overlayView.DrawRefColorDepiction(img)
    return overlayView.refmol2D, img

def getInput():
    fileChooser = JFileChooser()
    val = fileChooser.showSaveDialog(None)
    if val == JFileChooser.APPROVE_OPTION:
        f = fileChooser.selectedFile
        return f.getName()
    else:
        return None


class MyTranscoder(ImageTranscoder):
    def __init__(self):
        self.image = None

    def createImage(self, w, h):
        self.image = BufferedImage(w, h, BufferedImage.TYPE_INT_ARGB)
        return self.image

    def getImage(self):
        return self.image


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


class MolUtilities():
    def __init__(self):
        return

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

def getAstexColor(rowId):
    colors = ["gray","yellow","magenta","cyan","pink","orange"]
    if rowId <0:
        return colors[0]
    else:
        return colors[rowId%len(colors)]

class JobRunner():
    def __init__(self):
        self.progressMonitor = ProgressMonitor(vortex, "Progress", "Progress", 0, 100)
        self.progressMonitor.setMillisToPopup(0)
        self.progressMonitor.setMillisToDecideToPopup(0)
        self.isCancelled = False
        self.selectedColumn = None
        self.numHits = 100
        self.nameCB = None
        self.xmlrpc_client = XmlRpcClient()
        self.queryMol = None
        self.molUtilities = MolUtilities()
        return

    def show_interface(self):
        sketcher = MSketchPane()
        molUtilities = self.molUtilities
        if vtable is not None:
            selectedRows = vtable.getSelectedRows()
            if selectedRows.size() >= 1:
                selectedRow = selectedRows.get(0)
                mfm = vtable.getMolFileManager()
                queryMol1 = mfm.getMolFileAtRow(selectedRow)
                class MyRunnable(java.lang.Runnable):
                    def run(self):
                        mol = molUtilities.convertToChemaxonMol(queryMol1)
                        sketcher.setMol(mol)
                SwingUtilities.invokeLater(MyRunnable())

        sketcher.setPreferredSize(Dimension(800, 600))
        p = JPanel(BorderLayout())
        p.add(sketcher, BorderLayout.CENTER)
        p2 = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        self.dbCB = JComboBox(["Biogen In Stock (400K)","All Screening", "Charles River Diverse (660K)"])

        self.slider = JSlider(0,1000,50)
        self.slider.setMajorTickSpacing(200)
        self.slider.setPaintTicks(True)
        self.slider.setPaintLabels(True)
        self.slider.setLabelTable(self.slider.createStandardLabels(200))
        textField = JTextField(5)
        textField.setText("50")
        textField.setEditable(False)
        slider1 = self.slider
        class SliderChangeListener(ChangeListener):
            def stateChanged(self, e):
                textField.setText("%d"%slider1.getValue())

        self.slider.addChangeListener(SliderChangeListener())
        p2.add(JLabel("Database:"))
        p2.add(self.dbCB)
        p2.add(JLabel("No. of Solutions:"))
        p2.add(self.slider)
        p2.add(textField)
        p.add(p2,BorderLayout.SOUTH)
        p.setBorder(BorderFactory.createTitledBorder("Click to edit structure:"))

        ret = vortex.showInDialog(p, "Sketch input molecule for FastRocs Search")
        if ret == 0:
            if sketcher.getMol().isEmpty():
                return False
            else:
                self.queryMol = MolExporter.exportToFormat(sketcher.getMol(),"mol")
                print self.queryMol
                return True
        else:
            return False

    def run(self):
        this = self
        this.numHits = this.slider.getValue()
        if this.numHits == 0:
            this.numHits = 50

        xmlrpc_config = XmlRpcClientConfigImpl()

        if this.dbCB.getSelectedItem() == "Biogen In Stock (150K)":
            xmlrpc_config.setServerURL(URL(FASTROCS_SERVER_EVOTEC))
        elif this.dbCB.getSelectedItem() == "All Screening":
            xmlrpc_config.setServerURL(URL(FASTROCS_SERVER_ALL))
        elif this.dbCB.getSelectedItem() == "SoftFocus Diverse (660K)":
            xmlrpc_config.setServerURL(URL(FASTROCS_SERVER_SOFTFOCUSDIVERSE))
        else:
            vortex.alert("Unknown database!")
            return

        this.xmlrpc_client.setConfig(xmlrpc_config)

        class MyTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self, progressMonitor):
                self.progressMonitor = progressMonitor
                self.molFileManager = {}
                self.molUtilities = MolUtilities()
                SwingWorker.__init__(self)

            def doInBackground(self):
                self.super__publish([("Running FastROCS ...", 1)])
                molUtilities = self.molUtilities
                mol2d = molUtilities.convertToMol(this.queryMol)
                self.mol3d = molUtilities.get3DMol(mol2d)
                if self.mol3d is not None:
                    params = [String(self.mol3d).getBytes(), this.numHits, "sdf", "sdf"]
                    try:
                        qIdx = this.xmlrpc_client.execute("SubmitQuery", params)

                        while True:
                            statusParams = [qIdx, True]
                            result = this.xmlrpc_client.execute("QueryStatus", statusParams)
                            current = int(result[0])
                            total = int(result[1])
                            if total == 0:
                                continue
                            progress = 80 * current / total
                            self.super__publish([("Running ...", progress)])
                            if total <= current:
                                break

                        res = this.xmlrpc_client.execute("QueryResults", [qIdx])
                        res = StringUtil.fromBytes(res)
                        vortexTableModel = VortexTableModel()
                        structureColumn = vortexTableModel.findColumnWithName("MOLFILE", 1, Column.COLUMN_STRING)
                        imageColumn = vortexTableModel.findColumnWithName("Image", 1, Column.COLUMN_STRING)
                        nameColumn = vortexTableModel.findColumnWithName("Name", 1, Column.COLUMN_STRING)
                        shapeTanimotoColumn = vortexTableModel.findColumnWithName("ShapeTanimoto", 1,
                                                                                  Column.COLUMN_DOUBLE)
                        colorTanimotoColumn = vortexTableModel.findColumnWithName("ColorTanimoto", 1,
                                                                                  Column.COLUMN_DOUBLE)
                        comboTanimotoColumn = vortexTableModel.findColumnWithName("TanimotoCombo", 1,
                                                                                  Column.COLUMN_DOUBLE)
                        similarity2DColumn = vortexTableModel.findColumnWithName("Similairty 2D (MACCS)", 1,
                                                                                 Column.COLUMN_DOUBLE)
                        refMol3d = molUtilities.convertToMol(self.mol3d)
                        mol2d, images = getOverlayImagesRef(refMol3d)
                        structureColumn.setValueFromString(0, molUtilities.molToSDF(mol2d))
                        nameColumn.setValueFromString(0, "Reference molecule")
                        array = oedepict.oedepict.OEWriteImageToByteArray("svg", images)
                        imageColumn.setValueFromString(0, String(array))
                        shapeTanimotoColumn.setValueFromString(0, "1.0")
                        colorTanimotoColumn.setValueFromString(0, "1.0")
                        comboTanimotoColumn.setValueFromString(0, "2.0")
                        similarity2DColumn.setValueFromString(0, "1.0")
                        self.molFileManager[0] = self.mol3d
                        fingerprint_ref = oegraphsim.OEFingerPrint()
                        oegraphsim.oegraphsim.OEMakeMACCS166FP(fingerprint_ref, refMol3d)
                        #                molFileManager.addMoleFile(String(this.refSDF))
                        drugs = molUtilities.getMolListFromSDFString(res)
                        for id, oe_mol in enumerate(drugs):
                            print "%d out of %d" % (id, len(drugs))
                            print "processing mol :"
                            print  molUtilities.molToSDF(oe_mol)
                            if not id:
                                continue
                            self.super__publish([("building table ...", 80 + 20 * id / len(drugs))])
                            fingerprint_mol = oegraphsim.OEFingerPrint()
                            oegraphsim.oegraphsim.OEMakeMACCS166FP(fingerprint_mol, oe_mol)
                            similarity = oegraphsim.oegraphsim.OETanimoto(fingerprint_ref, fingerprint_mol)

                            self.molFileManager[id] = molUtilities.molToSDF(oe_mol)
                            #                    molFileManager.addMoleFile(String(sdf))
                            mol2d, overlay_images = getOverlayImages(refMol3d, oe_mol)
                            structureColumn.setValueFromString(id, molUtilities.molToSDF(mol2d))
                            nameColumn.setValueFromString(id, oe_mol.GetTitle())
                            array = oedepict.oedepict.OEWriteImageToByteArray("svg", overlay_images)
                            imageColumn.setValueFromString(id, String(array))
                            shapeTanimotoColumn.setValueFromString(id, oechem.oechem.OEGetSDData(oe_mol,
                                                                                                 shapeTanimotoColumn.getName()))
                            colorTanimotoColumn.setValueFromString(id, oechem.oechem.OEGetSDData(oe_mol,
                                                                                                 colorTanimotoColumn.getName()))
                            comboTanimotoColumn.setValueFromString(id, oechem.oechem.OEGetSDData(oe_mol,
                                                                                                 comboTanimotoColumn.getName()))
                            similarity2DColumn.setValueFromString(id, "%5.2f" % similarity)

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
                    viewPanel = JPanel(BorderLayout())
                    viewPanel.add(molViewer, BorderLayout.CENTER)

                    scriptField = JTextField(50)
                    scriptBtn = JButton("Run Command:")

                    class ViewerBtnActionListener(ActionListener):
                        def actionPerformed(self, e):
                            molViewer.execute(scriptField.getText())
                            molViewer.dirtyRepaint()

                    scriptBtn.addActionListener(ViewerBtnActionListener())

                    class ViewerKeyListener(KeyListener):
                        def keyTyped(self, e):
                            dialog = viewPanel.getParent()
                            if dialog is not None:
                                dialog.getRootPane().setDefaultButton(scriptBtn)
                            return

                        def keyPressed(self, e):
                            if e.getKeyCode() == KeyEvent.VK_ENTER:
                                scriptBtn.doClick()
                            return

                        def keyReleased(self, e):
                            return

                    scriptField.addKeyListener(ViewerKeyListener())

                    viewerBtnPanel = JPanel(FlowLayout(FlowLayout.CENTER, 20, 10))
                    viewerBtnPanel.add(scriptBtn)
                    viewerBtnPanel.add(scriptField)
                    #                    viewPanel.add(viewerBtnPanel,BorderLayout.PAGE_END)

                    visibleMolDict = {}
                    molDict = {}
                    molFileManager = self.molFileManager
                    molUtilities = self.molUtilities

                    class FastRocsTable(VortexJTable2):
                        def vortexTableChanged(self, tableModel, action, i1, i2, o):
                            if action == "SelectionChanged" and viewPanel.isShowing():
                                rowids = vortexTableModel.getSelectedRows().getArray()
                                for rowId in rowids:
                                    if not visibleMolDict.has_key(rowId):
                                        if not molDict.has_key(rowId):
                                            molFile = String(molFileManager[rowId])
                                            oemol = molUtilities.convertToMol(molFile)
                                            oechem.oechem.OESuppressHydrogens(oemol, True)
                                            molString = AstexFile(
                                                ByteArrayInputStream(String(molUtilities.molToSDF(oemol)).getBytes()))
                                            molecule = MoleculeIO.readMDLMol(molString)
                                            molDict[rowId] = molecule
                                            molecule.setName("molecule_at_row_%d" % rowId)
                                            molViewer.addMolecule(molecule)
                                        else:
                                            molecule = molDict[rowId]

                                        molecule.setDisplayed(1)
                                        visibleMolDict[rowId] = True
                                        molViewer.execute("display lines off molecule %s;"%molecule.getName())
                                        molViewer.execute("display cylinders on molecule %s;"%molecule.getName())
                                        molViewer.execute("color %s molecule %s and element 6;"%(getAstexColor(rowId),molecule.getName()))

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
                                            molDict[rowId].setDisplayed(0)

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
                    # btn2 = JButton("View Selected In PyMol")
                    # btnPanel.add(btn2)
                    # multiplexCB = JCheckBox("Single Obj?",False)
                    # btnPanel.add(multiplexCB)
                    #
                    # btn3 = JButton("Reset Pymol Path")
                    # btnPanel.add(btn3)
                    btn4 = JButton("Export to SDF")
                    btnPanel.add(btn4)

                    mainPanel.add(btnPanel, BorderLayout.PAGE_END)
                    ws = vortex.createWorkspace(vortexTableModel,
                                                "%s_%s" % ("Fast Rocs Result", vortex.getWorkspaceCount()))
                    ws.addComponent(mainPanel, "FastRocs")
                    ws.updateLayout(VortexWorkspace.LAYOUT_GRID)

                    class View3DActionListener(ActionListener):
                        def actionPerformed(self, e):
                            if not viewPanel.isShowing():
                                ws.addComponent(viewPanel, "3D Viewer")
                                ws.updateLayout(VortexWorkspace.LAYOUT_GRID)
                                table.vortexTableChanged(vortexTableModel, "SelectionChanged", 0, 0, None)
                                #                            vortex.showInDialog(viewPanel,"molecules")

                    btn.addActionListener(View3DActionListener())

                    class ViewPyMolActionListener(ActionListener):
                        def actionPerformed(self, e):
                            pymolClient = PymolClient()
                            if pymolClient.isSuccess():
                                tmpFile = Vortex.getTempFile("sdf")
                                f = open(tmpFile, "w")
                                rowids = vortexTableModel.getSelectedRows().getArray()
                                for rowId in rowids:
                                    molFile = molFileManager[rowId]
                                    f.write(molFile)
                                f.close()
                                if multiplexCB.isSelected():
                                    multiplex = 0
                                else:
                                    multiplex = 1
                                pymolClient.open_file(tmpFile, multiplex)

                    # btn2.addActionListener(ViewPyMolActionListener())

                    class ResetPymolPathListener(ActionListener):
                        def actionPerformed(self, e):
                            vortex.deleteProperty("pymolexe")
                            vortex.deleteProperty("pymolpath")
                            vortex.alert("Pymol location is reset.")

                    # btn3.addActionListener(ResetPymolPathListener())


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
                except ExecutionException, e:
                    traceback.print_exc()
                finally:
                    self.progressMonitor.close()
                return

        task = MyTask(self.progressMonitor)
        task.execute()


jobRunner = JobRunner()
if jobRunner.show_interface():
    jobRunner.run()
