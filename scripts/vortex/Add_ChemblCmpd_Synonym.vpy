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
import javax.swing.JTabbedPane as JTabbedPane
import javax.swing.JTextArea as JTextArea
import javax.swing.event.ChangeListener as ChangeListener
import java.lang.Float as Float
import java.util.concurrent.ExecutionException as ExecutionException
import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.table.VortexJTable2 as VortexJTable2
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.VortexApplication as VortexApplication
from com.dotmatics.chemlink import ChemLink


sys.path.append(os.path.join(vortex.getVortexFolder(),"openeye"))

from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

import openeye.oechem as oechem
import uk.ac.ebi.chemblws as chemblws
import uk.ac.ebi.chemblws.restclient.ChemblRestClient as ChemblRestClient
import uk.ac.ebi.chemblws.domain.Target
import org.springframework.context.support.ClassPathXmlApplicationContext as ApplicationContext

class MolUtilities():
    def __init__(self):
        return

    def convertToMol(self, molString, format):
        ifs = oechem.oemolistream()
        ifs.SetFormat(format)
        ifs.openstring(molString)
        mol = oechem.OEGraphMol()
        oechem.oechem.OEReadMolecule(ifs,mol)
        return mol

    def convertToInchiKey(self,mol):
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat.INCHIKEY)
        ofs.openstring()
        oechem.oechem.OEWriteMolecule(ofs,mol)
        return ofs.GetString()

    def molToSmiles(self, mol):
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat.SMI)
        ofs.openstring()
        oechem.oechem.OEWriteMolecule(ofs,mol)
        return ofs.GetString()

class JobRunner():
    def __init__(self):
        self.progressMonitor = ProgressMonitor(vortex,"Progress","Progress", 0, 100)
        self.progressMonitor.setMillisToPopup(0)
        self.progressMonitor.setMillisToDecideToPopup(0)
        self.isCancelled = False
        self.selectedColumn = None
        self.numHits = 100
        self.nameCB = None
        ##################
        # from java.lang import ClassLoader
        # cl = ClassLoader.getSystemClassLoader()
        # paths = map(lambda url: url.getFile(), cl.getURLs())
        # file_io = open("/Users/jfeng1/junk.log", "w")
        # for p in paths:
        #     file_io.write(p+"\n")
        # file_io.close()
        ##################
        applicationContext = ApplicationContext("applicationContext.xml")
        self.chemblClient = applicationContext.getBean("chemblRestClient")
        return

    def show_interface(self):
        structureCols = []
        for colName in vtable.getColumnNames(0): #not showing type
            col = vtable.findColumnWithName(colName,0)
            if col is not None and col.isStructureColumn():
                structureCols.append(colName)
        if len(structureCols) > 0:
            content = javax.swing.JPanel()
            self.structureCB = JComboBox(structureCols)
            p = JPanel()
            p.add(JLabel("Select Structure Column:"))
            p.add(self.structureCB)
            layout.fill(content,p,0,0)
            ret = not vortex.showInDialog(content, "Pick structure column:")
            if ret:
                self.isCancelled = False
            else:
                self.isCancelled = True
                return False

            colName = self.structureCB.getSelectedItem()
            self.selectedColumn = vtable.findColumnWithName(colName,0)
            if self.selectedColumn is not None:
                return True
            else:
                vortex.alert("Column is not available.")
                return False
        else:
            vortex.alert("No structural columns available.")
            self.isCancelled = True
            return False

    def run(self):
        this = self

        class MyTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self, progressMonitor):
                self.progressMonitor = progressMonitor
                SwingWorker.__init__(self)

            def doInBackground(self):
                topHits = []
                mfm = vtable.getMolFileManager()
                nRows = vtable.getRealRowCount()
                rows = range(0,nRows)

                molUtilities = MolUtilities()
                chemblCol = vtable.findColumnWithName("Compound Name",1)
                for row in rows:
                    message = "Processing mol %d..."%row
                    progress = 100*row/nRows
                    self.super__publish([(message, progress)])
                    mol = None
                    molString = None
                    colName = this.selectedColumn.getName()
                    if colName =="Structure":
                        molString = mfm.getMolFileAtRow(row)
                        mol = molUtilities.convertToMol(molString,oechem.OEFormat.MDL)
                    else:
                        molString = this.selectedColumn.getValueAsString(row)
                        mol = molUtilities.convertToMol(molString,oechem.OEFormat.SMI)
                    inchiKey = molUtilities.convertToInchiKey(mol)
                    print inchiKey
                    try:
                        compound = this.chemblClient.getCompoundByStdInChiKey(inchiKey)
                        # import uk.ac.ebi.chemblws.domain.Compound as Compound
                        # compound = Compound()
                        compound_name = None
                        if compound.getPreferredCompoundName() != None and len(compound.getPreferredCompoundName())>0:
                            compound_name = compound.getPreferredCompoundName()
                        elif compound.getSynonyms() != None and len(compound.getSynonyms())>0:
                            compound_name = compound.getPreferredCompoundName()
                        else:
                            compound_name = ""

                        chemblCol.setValueFromString(row,compound_name)
                    except:
                        traceback.print_exc()
                chemblCol.calculateWidth()
                return None

            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setNote(pair[0])
                self.progressMonitor.setProgress(pair[1])

            def done(self):
                try:
                    self.get()  #raise exception if abnormal completion
                    vtable.fireTableStructureChanged()
                except ExecutionException, e:
                    traceback.print_exc()
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
                break
