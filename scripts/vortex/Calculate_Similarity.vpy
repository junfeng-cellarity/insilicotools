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
import openeye.oegraphsim as oegraphsim

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
        self.isCancelled = True
        self.selectedColumn = None
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
        this = self
        chemLink = ChemLink()

        if vtable is not None:
            selectedRows = vtable.getSelectedRows()

            if selectedRows.size() >= 1:
                selectedRow = selectedRows.get(0)
                mfm = vtable.getMolFileManager()
                queryMol1 = mfm.getMolFileAtRow(selectedRow)
                class MyRunnable(java.lang.Runnable):
                    def run(self):
                        chemLink.setMolFile(queryMol1)
                        chemLink.repaint()

                SwingUtilities.invokeLater(MyRunnable())
            chemLink.setPreferredSize(Dimension(600, 400))

        p = JPanel(BorderLayout())
        p.setBorder(BorderFactory.createTitledBorder("Click to edit structure:"))
        p.add(chemLink,BorderLayout.CENTER)

        typeCB = JComboBox(["MACCS166","Circular","Lingo"])
        p2 = JPanel()
        p2.add(JLabel("Simiarity Type:"))
        p2.add(typeCB)
        p.add(p2,BorderLayout.SOUTH)

        dialog = JDialog(VortexApplication.vortexApplication,"Calculate Similarity",True)
        dialog.setSize(Dimension(800,600))
        dialog.setLocationRelativeTo(None)
        dialog.setDefaultCloseOperation(JDialog.HIDE_ON_CLOSE)

        btnPanel = JPanel()
        okButton = JButton("Ok")
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

        btnPanel.add(okButton)
        btnPanel.add(cancelButton)

        mainPanel = JPanel(BorderLayout())
        mainPanel.add(p,BorderLayout.CENTER)
        mainPanel.add(btnPanel,BorderLayout.SOUTH)
        dialog.getContentPane().add(mainPanel)
        dialog.setVisible(True)
        if typeCB.getSelectedItem() == "MACCS166":
            self.fingerprintType = oegraphsim.OEFPType.MACCS166
        elif typeCB.getSelectedItem() == "Circular":
            self.fingerprintType = oegraphsim.OEFPType.Circular
        else:
            self.fingerprintType = oegraphsim.OEFPType.Lingo

        if not this.isCancelled:
            molfile = str(chemLink.getMolFile())
            if len(molfile.strip()) == 0:
                return False
            else:
                molUtilities = MolUtilities()
                self.inputMol = molUtilities.convertToMol(molfile,oechem.OEFormat.MDL)
                if self.inputMol.NumAtoms()>0:
                    return True
                else:
                    return False
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

                fp = oegraphsim.OEFingerPrint()
                oegraphsim.oegraphsim.OEMakeFP(fp, this.inputMol, this.fingerprintType)
                simColumn  = vtable.findColumnWithName("SimilarityScore",1,Column.COLUMN_DOUBLE)
                for row in rows:
                    message = "Processing mol %d..."%(row+1)
                    progress = 100*row/nRows
                    self.super__publish([(message, progress)])
                    mfm = vtable.getMolFileManager()
                    molFile = mfm.getMolFileAtRow(row)
                    if molFile is None or len(molFile) == 0:
                        continue
                    mol = molUtilities.convertToMol(molFile,oechem.OEFormat.MDL)
                    fp2 = oegraphsim.OEFingerPrint()
                    oegraphsim.oegraphsim.OEMakeFP(fp2, mol, this.fingerprintType)
                    tanimoto = oegraphsim.oegraphsim.OETanimoto(fp,fp2)
                    simColumn.setValueFromString(row,"%5.2f"%tanimoto)

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
                    vortex.alert("Failed to calculate similarity.")
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
