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
import java.util.Random as Random
import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.table.VortexJTable2 as VortexJTable2
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.VortexApplication as VortexApplication



sys.path.append(os.path.join(vortex.getVortexFolder(),"openeye"))

from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

import openeye.oechem as oechem
import openeye.oegraphsim as oegraphsim
import org.apache.commons.math3.stat.clustering.KMeansPlusPlusClusterer as KMeansPlusPlusClusterer
import org.apache.commons.math3.stat.clustering.Clusterable as Clusterable

class Compound(Clusterable):
    def __init__(self, id, fp):
        self.id = id
        self.fp = fp
        pass

    def get_distance(self,a,b):
        return 1.0 - oegraphsim.oegraphsim.OETanimoto(a.fp,b.fp)

    def distanceFrom(self,a):
        distance = 1.0 - oegraphsim.oegraphsim.OETanimoto(self.fp,a.fp)
        return distance

    def centroidOf(self, group):
        if group is None or len(group) == 0:
            return self
        if len(group):
            return group[0]

        name_dict = {}
        for a in group:
            name_dict[a.id] = a

        distance_dict = {}
        n = len(group) - 1
        for a in group:
            sum = 0
            for b in group:
                if a.id != b.id:
                    sum += self.get_distance(a,b)
            distance_dict[a.id] = sum/n
        return sorted(group,key=lambda cmpd:distance_dict[cmpd.id])[0]

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
        self.numClusters = 0
        self.fingerprintType = None
        self.clusterMethod = None
        self.fingerprintTypeName = None
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
        selectedRows = None
        if vtable is not None:
            selectedRows = vtable.getSelectedRows()
            if selectedRows.size()<5:
                vortex.alert("no or too few (<5) molecules selected.")
                self.isCancelled = True
                return False
        else:
            vortex.alert("No Workspace available.")
            self.isCancelled = True
            return False

        typeCB = JComboBox(["MACCS166","Circular","Lingo","Path","Tree"])
        p1 = JPanel()
        p1.add(JLabel("Fingerprint:"))
        p1.add(typeCB)

        p2 = JPanel()
        p2.add(JLabel("Number of Clusters:"))
        numClustersField = JTextField(10)
        numClusters = int(selectedRows.size())/10
        if numClusters <2:
            numClusters = 2
        numClustersField.setText("%d"%numClusters)
        p2.add(numClustersField)

        clusterCB = JComboBox(["KMean"])
        p3 = JPanel()
        p3.add(JLabel("Clustering Method:"))
        p3.add(clusterCB)

        numTrialsCB = JComboBox(['5','10','20','50','100'])
        p4 = JPanel()
        p4.add(JLabel("Num. of Trials:"))
        p4.add(numTrialsCB)

        p0 = JPanel(GridLayout(4,1))
        p0.add(p1)
        p0.add(p2)
        p0.add(p3)
        p0.add(p4)



        dialog = JDialog(VortexApplication.vortexApplication,"Clustering",True)
        dialog.setSize(Dimension(300,250))
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
        mainPanel.add(p0,BorderLayout.CENTER)
        mainPanel.add(btnPanel,BorderLayout.SOUTH)
        dialog.getContentPane().add(mainPanel)
        dialog.setVisible(True)


        if not this.isCancelled:
            try:
                self.numClusters = int(numClustersField.getText())
            except:
                vortex.alert("Invalid cluster number.")
                self.isCancelled = False
                return False

            self.numTrials = int(numTrialsCB.getSelectedItem())
            self.clusterMethod = clusterCB.getSelectedItem()

            if typeCB.getSelectedItem() == "MACCS166":
                self.fingerprintType = oegraphsim.OEFPType.MACCS166
            elif typeCB.getSelectedItem() == "Circular":
                self.fingerprintType = oegraphsim.OEFPType.Circular
            elif typeCB.getSelectedItem() == "Path":
                self.fingerprintType = oegraphsim.OEFPType.Path
            elif typeCB.getSelectedItem() == "Tree":
                self.fingerprintType = oegraphsim.OEFPType.Tree
            elif typeCB.getSelectedItem() == "Lingo":
                self.fingerprintType = oegraphsim.OEFPType.Lingo
            else:
                vortex.alert("Unknown fingerprint type!")
                self.isCancelled = False
                return False
            self.fingerprintTypeName = typeCB.getSelectedItem()
        else:
            return False

        return True

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
                selectedRows = vtable.getSelectedRows()
                nRows = selectedRows.size()
                molUtilities = MolUtilities()

                mfm = vtable.getMolFileManager()
                clusterableList = []
                for row in selectedRows.getArray():
                    message = "Processing mol %d..."%(row+1)
                    progress = 50*row/nRows
                    self.super__publish([(message, progress)])
                    molFile = mfm.getMolFileAtRow(row)
                    if molFile is None or len(molFile) == 0:
                        continue
                    mol = molUtilities.convertToMol(molFile,oechem.OEFormat.MDL)
                    if mol is None or mol.NumAtoms()==0:
                        continue
                    fp = oegraphsim.OEFingerPrint()
                    oegraphsim.oegraphsim.OEMakeFP(fp, mol, this.fingerprintType)
                    clusterableList.append(Compound(row,fp))
                clusterer = KMeansPlusPlusClusterer(Random())
                clusters = clusterer.cluster(clusterableList,this.numClusters,this.numTrials,-1)
                clusters = sorted(clusters, key=lambda cluster:cluster.getPoints().size(), reverse=True)
                clusterColumn = vtable.findColumnWithName("KMeanCluster_%s"%this.fingerprintTypeName,1, Column.COLUMN_INT)
                clusterCenterColumn = vtable.findColumnWithName("KMeanClusterCenter_%s"%this.fingerprintTypeName,1, Column.COLUMN_INT)
                for clusterId,cluster in enumerate(clusters):
                    clusterCenterColumn.setValueFromString(cluster.getCenter().id,"1")
                    for c in cluster.getPoints():
                        clusterColumn.setValueFromString(c.id,"%d"%clusterId)
                        if c.id != cluster.getCenter().id:
                            clusterCenterColumn.setValueFromString(c.id,"0")
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
                    vortex.alert("Failed to run clustering.")
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
