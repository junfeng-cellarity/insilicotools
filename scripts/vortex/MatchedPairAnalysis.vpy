import sys
import math
import traceback
import os
import javax.swing.ProgressMonitor as ProgressMonitor
import javax.swing.SwingWorker as SwingWorker
import javax.swing.JComboBox as JComboBox
import javax.swing.JTextField as JTextField
import javax.swing.JButton as JButton
import javax.swing.JDialog as JDialog
import java.awt.event.ActionListener as ActionListener
import javax.swing.event.DocumentListener as DocumentListener
import javax.swing.JLabel as JLabel
import javax.swing.JPanel as JPanel
import javax.swing.SwingUtilities as SwingUtilities
import javax.swing.BorderFactory as BorderFactory
import java.awt.BorderLayout as BorderLayout
import java.awt.GridLayout as GridLayout
import java.awt.Dimension as Dimension
import javax.swing.JScrollPane as JScrollPane
import javax.swing.JSplitPane as JSplitPane
import javax.swing.JPanel as JPanel
import java.util.concurrent.ExecutionException as ExecutionException
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
import java.lang.StringBuffer as StringBuffer

import java.beans.PropertyChangeListener as PropertyChangeListener

import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.table.VortexJTable2 as VortexJTable2
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.VortexApplication as VortexApplication
import com.dotmatics.vortex.plots.VortexScatterPlot as vp

from com.dotmatics.chemlink import ChemLink

sys.path.append(os.path.join(vortex.getVortexFolder(),"openeye"))

from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

import openeye.oeomega as oeomega
import openeye.oechem as oechem
import openeye.oegraphsim as oegraphsim
import openeye.oechem.OEGraphMol as OEGraphMol
import chemaxon.marvin.beans.MSketchPane as MSketchPane
import chemaxon.formats.MolImporter as MolImporter
import chemaxon.formats.MolExporter as MolExporter
import chemaxon.struc.RxnMolecule as RxnMolecule
import com.chemaxon.mapper.AutoMapper as AutoMapper
import com.dotmatics.vortex.plots.VortexScatterPlot as ScatterPlot

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

    def molToSmiles(self, mol, stereo):
        try:
            if stereo:
                return oechem.oechem.OEMolToSmiles(mol).split()[0]
            else:
                smiles = StringBuffer()
                oechem.oechem.OECreateCanSmiString(smiles,mol)
                return smiles.toString().split()[0]
        except:
            return None

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

class Enumerator():
    def __init__(self,smirks):
        self.libraryGen = oechem.OELibraryGen()
        self.libraryGen.Init(smirks)
        self.libraryGen.SetExplicitHydrogens(False)
        self.libraryGen.SetValenceCorrection(True)
        self.libraryGen.SetRemoveUnmappedFragments(True)

    def addReagent(self,mol):
        oechem.oechem.OEAssignAromaticFlags(mol)
        self.libraryGen.AddStartingMaterial(mol,0)

    def enumerate(self):
        molList = []
        if self.libraryGen.IsValid():
            productIter = self.libraryGen.GetProducts()
            while productIter.hasNext():
                product = productIter.next()
                if product is not None and product.NumAtoms() >0:
                    molList.append(OEGraphMol(product))
        return molList

class JobRunner():
    def __init__(self):
        self.progressMonitor = ProgressMonitor(vortex,"Progress","Progress", 0, 100)
        self.progressMonitor.setMillisToPopup(0)
        self.progressMonitor.setMillisToDecideToPopup(0)
        self.isCancelled = True
        self.selectedColumn = None
        self.molUtilities = MolUtilities()
        ##################
        # from java.lang import ClassLoader
        # cl = ClassLoader.getSystemClassLoader()
        # paths = map(lambda url: url.getFile(), cl.getURLs())
        # file_io = open("/Users/jfeng1/junk.log", "w")
        # for p in paths:
        #     file_io.write(p+"\n")
        # file_io.close()
        ##################
        return

    def show_interface(self):
        sketcher = MSketchPane()
        molUtilities = self.molUtilities
        autoMapper = AutoMapper()
        autoMapper.setMappingStyle(AutoMapper.MappingStyle.MATCHING)
        self.columns = []
        if vtable is not None:
            mmp_mol = vortex.getProperty("mmp_mol")
            if mmp_mol is not None and len(mmp_mol)>0:
                sketcher.setMol(MolImporter.importMol(mmp_mol))
            else:
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


            for column in vtable.getColumnsArray():
                if column.getForceNumeric() or column.isNumeric():
                    self.columns.append(column.getName())

        self.columnCB = JComboBox(self.columns)
        self.chiralCB = JCheckBox("Chiral Sensitive?")

        sketcher.setPreferredSize(Dimension(800, 600))
        p = JPanel(BorderLayout())
        p.add(sketcher, BorderLayout.CENTER)
        p2 = JPanel(FlowLayout(FlowLayout.CENTER,20,10))

        p2.add(self.chiralCB)
        p2.add(JLabel("Data Colomn:"))
        p2.add(self.columnCB)

        autoMapBtn = JButton("Auto Map")
        class AutoMapListener(ActionListener):
            def actionPerformed(self,e):
                mol = sketcher.getMol()
                rxnMol = RxnMolecule.getReaction(mol)
                if not rxnMol.isEmpty():
                    autoMapper.map(rxnMol)
                    sketcher.setMol(rxnMol)
        autoMapBtn.addActionListener(AutoMapListener())
        p2.add(autoMapBtn)

        p.add(p2,BorderLayout.SOUTH)
        p.setBorder(BorderFactory.createTitledBorder("Use Reaction Arrow To Draw Reaction"))

        ret = vortex.showInDialog(p, "Sketch transforming reaciton")
        if ret == 0:
            if sketcher.getMol().isEmpty():
                return False
            else:
                self.smirks = sketcher.getMol("smarts:ah")
                vortex.setProperty("mmp_mol",sketcher.getMol("mol"))
                mol = sketcher.getMol()
                rxnMol = RxnMolecule.getReaction(mol)
                self.reactantSmi = rxnMol.getReactant(0).exportToFormat("smiles:a")
                self.productSmi = rxnMol.getProduct(0).exportToFormat("smiles:a")
                return True
        else:
            return False


    def run(self):

        this = self
        class MyTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self, progressMonitor):
                self.progressMonitor = progressMonitor
                self.progressMonitor.setMillisToDecideToPopup(0)
                self.progressMonitor.setMillisToPopup(0)
                self.progressMonitor.setNote("Calculating similarity ...")
                self.progressMonitor.setProgress(1)
                SwingWorker.__init__(self)

            def doInBackground(self):
                molUtilities = MolUtilities()

                nRows = vtable.getRealRowCount()
                rows = range(0,nRows)
                mol_dict ={}
                pair_dict = {}
                product_dict = {}
                n = 0
                for column in vtable.getColumnsArray():
                    if column.getName().startswith("MMP Smiles"):
                        n += 1

                # mmpRoleColumn  = vtable.findColumnWithName("MMP Smiles %d"%n)
                # if mmpRoleColumn is not None:
                #     vtable.deleteColumn(mmpRoleColumn)
                # mmpSeqColumn = vtable.findColumnWithName("MMP Id %d"%n)
                # if mmpSeqColumn is not None:
                #     vtable.deleteColumn(mmpSeqColumn)

                self.mmp_smiles_column_name = "MMP Smiles %d" % n
                mmpRoleColumn = vtable.findColumnWithName(self.mmp_smiles_column_name, True, Column.COLUMN_STRING)
                #mmpRoleColumn  = vtable.findColumnWithName(self.mmp_smiles_column_name,False)
                self.mmp_seq_column_name = "MMP Id %d" % n
                mmpSeqColumn = vtable.findColumnWithName(self.mmp_seq_column_name, True, Column.COLUMN_INT)
                #mmpSeqColumn = vtable.findColumnWithName(self.mmp_seq_column_name,False)
                if mmpRoleColumn and mmpSeqColumn:
                    for row in rows:
                        mmpRoleColumn.setValueFromString(row,None)
                        mmpSeqColumn.setValueFromString(row,None)
                        message = "Processing mol %d..."%(row+1)
                        progress = 100*row/nRows
                        self.super__publish([(message, progress)])
                        mfm = vtable.getMolFileManager()
                        molFile = mfm.getMolFileAtRow(row)
                        if molFile is None or len(molFile) == 0:
                            continue
                        mol = molUtilities.convertToMol(molFile)
                        uniq_smi = molUtilities.molToSmiles(mol,this.chiralCB.isSelected())
                        if uniq_smi is None:
                            continue
                        enumerator = Enumerator(this.smirks)
                        enumerator.addReagent(mol)
                        products = enumerator.enumerate()
                        mol_dict[uniq_smi] = row
                        if len(products)>0:
                            pair_dict[uniq_smi] = 1
                            product_smi = molUtilities.molToSmiles(products[0],this.chiralCB.isSelected())
                            if product_smi is not None:
                                product_dict[product_smi] = row

                    if product_dict is not None and len(product_dict) >0:
                        mmp_rows = []
                        for product_smi in product_dict.keys():
                            row = product_dict[product_smi]
                            if product_smi in mol_dict.keys():
                                mmp_rows.append(row)

                        mmp_rows.sort()
                        row_dict = {}
                        for id,row in enumerate(mmp_rows):
                            row_dict[row] = id

                        for product_smi in product_dict.keys():
                            if product_smi in mol_dict.keys():
                                row = product_dict[product_smi]
                                mmpSeqColumn.setValueFromString(row, str(row_dict[row]))
                                mmpRoleColumn.setValueFromString(row,this.reactantSmi)
                                mmpSeqColumn.setValueFromString(mol_dict[product_smi],str(row_dict[row]))
                                mmpRoleColumn.setValueFromString(mol_dict[product_smi],this.productSmi)

            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setNote(pair[0])
                self.progressMonitor.setProgress(pair[1])

            def done(self):
                try:
                    self.get()  #raise exception if abnormal completion
                    mmpRoleColumn = vtable.findColumnWithName(self.mmp_smiles_column_name,False)
                    mmpSeqColumn = vtable.findColumnWithName(self.mmp_seq_column_name,False)
                    dataColumnName = this.columnCB.getSelectedItem()
                    dataColumn = vtable.findColumnWithName(dataColumnName, False)
                    plt = None
                    for k in range(vws.plotCount):
                        if isinstance(vws.getPlot(k),ScatterPlot):
                            plt = vws.getPlot(k)
                            break
                    if plt:
                        plt.setXColumn(vtable.getColumnId(mmpRoleColumn))
                        plt.setYColumn(vtable.getColumnId(dataColumn))
                        plt.setColourColumn(vtable.getColumnId(mmpRoleColumn))
                        plt.setSymbolColumn(vtable.getColumnId(mmpRoleColumn))
                        plt.setDisplayLines(True)
                        plt.setLineGroupColumn(vtable.getColumnId(mmpSeqColumn))
                        plt.setShowStructuresOnAxes(True)
                        vtable.fireTableStructureChanged()

                except ExecutionException, e:
                    traceback.print_exc()
                    vortex.alert(traceback.format_exc())
                finally:
                    self.progressMonitor.close()
                return

        task = MyTask(self.progressMonitor)
        task.execute()

if vws is None:
    vortex.alert("You must have a workspace loaded...")
else:
    jobRunner = JobRunner()
    while 1:
        if jobRunner.show_interface():
            jobRunner.run()
            break
        else:
            if jobRunner.isCancelled:
                vortex.alert("User cancelled.")
                break
            else:
                vortex.alert("No molecule available.")
