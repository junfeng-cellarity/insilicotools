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

descriptor_dict = {
    "Circular":oegraphsim.OEFPType.Circular,
    "Lingo":oegraphsim.OEFPType.Lingo,
    "MACSS":oegraphsim.OEFPType.MACCS166,
    "Path":oegraphsim.OEFPType.Path,
    "Tree":oegraphsim.OEFPType.Tree,
}

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

        typeCB = JComboBox(descriptor_dict.keys())
        p = JPanel()
        p.add(JLabel("Simiarity Type:"))
        p.add(typeCB)

        dialog = JDialog(VortexApplication.vortexApplication,"Select descriptor:",True)
        dialog.setSize(Dimension(300,100))
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
        self.fingerprintType = descriptor_dict.get(typeCB.getSelectedItem())

        if not this.isCancelled:
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
                self.progressMonitor.setNote("Generating fingerprint ...")
                self.progressMonitor.setProgress(1)
                SwingWorker.__init__(self)

            def doInBackground(self):
                molUtilities = MolUtilities()

                nRows = vtable.getRealRowCount()
                rows = range(0,nRows)

                fp_columns = []
                for row in rows:
                    message = "Processing mol %d..."%(row+1)
                    progress = 100*row/nRows
                    self.super__publish([(message, progress)])
                    mfm = vtable.getMolFileManager()
                    molFile = mfm.getMolFileAtRow(row)
                    if molFile is None or len(molFile) == 0:
                        continue
                    mol = molUtilities.convertToMol(molFile,oechem.OEFormat.MDL)
                    fp = oegraphsim.OEFingerPrint()
                    oegraphsim.oegraphsim.OEMakeFP(fp, mol, this.fingerprintType)
                    if row ==0:
                        for colIdx in range(0,fp.GetSize()):
                            col = vtable.findColumnWithName("%d"%colIdx, 1, Column.COLUMN_INT)
                            fp_columns.append(col)
                    for idx,col in enumerate(fp_columns):
                        if fp.IsBitOn(idx):
                            col.setValueFromString(row,"1")
                        else:
                            col.setValueFromString(row,"0")

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
