import sys
import math
import traceback
import os
import javax.swing.ProgressMonitor as ProgressMonitor
import javax.swing.SwingWorker as SwingWorker
import javax.swing.JComboBox as JComboBox
import javax.swing.JLabel as JLabel
import javax.swing.JPanel as JPanel
import javax.swing.JScrollPane as JScrollPane
import java.util.concurrent.ExecutionException as ExecutionException
import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.table.VortexJTable2 as VortexJTable2
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace

sys.path.append(os.path.join(vortex.getVortexFolder(),"openeye"))

from initOEChem import VortexInit
VortexInit(vortex).initOpeneye()
import openeye.oechem as oechem
import openeye.oegraphsim as oegraphsim

class MolUtilities():
    def __init__(self):
        return

    def convertToMol(self, molString):
        ifs = oechem.oemolistream()
        ifs.SetFormat(oechem.OEFormat.MDL)
        ifs.openstring(molString)
        mol = oechem.OEGraphMol()
        oechem.oechem.OEReadMolecule(ifs,mol)
        return mol

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
        self.numHits = 1000
        self.nameCB = None
        return

    def show_interface(self):
        cols = []
        nameCols = []
        for colName in vtable.getColumnNames(0): #not showing type
            col = vtable.findColumnWithName(colName,0)
            if col.getForceNumeric() or col.isNumeric():
                cols.append(colName)
            else:
                if not col.isStructureColumn():
                    nameCols.append(colName)
        if len(cols) > 0:
            content = javax.swing.JPanel()
            self.duallist = components.VortexDualList()
            self.duallist.setItems(cols, None)
            layout.fill(content, self.duallist, 0, 0)
            if len(nameCols) > 0:
                self.nameCB = JComboBox(nameCols)
                p = JPanel()
                p.add(JLabel("Select Name Column:"))
                p.add(self.nameCB)
                layout.fill(content,p,0,1)
            ret = not vortex.showInDialog(content, "Choose one column ONLY:")
            if ret:
                self.isCancelled = False
            else:
                self.isCancelled = True
                return False

            if len(self.duallist.getItems()) == 1:
                colIdx = self.duallist.getItems()[0]
                colName = cols[colIdx]
                self.selectedColumn = vtable.findColumnWithName(colName,0)
                if self.selectedColumn is not None:
                    return True
                else:
                    vortex.alert("Column is not available.")
                    return False
            else:
                if(len(self.duallist.getItems()) == 0):
                    vortex.alert("Need one column selected.")
                else:
                    vortex.alert("Only one column allowed!")
                return False
        else:
            vortex.alert("No column available")
            self.isCancelled = True
            return False

    def run(self):
        this = self
        nameCol = None
        if self.nameCB is not None:
            nameCol = vtable.findColumnWithName(self.nameCB.getSelectedItem(),0)

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

                molFPs = []
                smilesList = []
                molUtilities = MolUtilities()
                for row in rows:
                    molString = mfm.getMolFileAtRow(row)
                    mol = molUtilities.convertToMol(molString)
                    smiles = molUtilities.molToSmiles(mol)
                    fp = oegraphsim.OEFingerPrint()
                    oegraphsim.oegraphsim.OEMakeFP(fp, mol, oegraphsim.OEFPType.MACCS166)
                    molFPs.append(fp)
                    smilesList.append(smiles)

                for row1 in range(0,nRows-1):
                    message = "Processing mol %d..."%row1
                    progress = 100*row1/nRows
                    print "Progress %d"%progress
                    self.super__publish([(message, progress)])
                    strObj1 = this.selectedColumn.getValueAsString(row1)
                    name1 = "mol_%d"%row1
                    if nameCol is not None:
                        name1 = nameCol.getValueAsString(row1)
                    if strObj1 is not None:
                        try:
                            act1 = float(strObj1)
                            fp1 = molFPs[row1]
                            for row2 in range(row1+1,nRows):
                                fp2 = molFPs[row2]
                                strObj2 = this.selectedColumn.getValueAsString(row2)
                                name2 = "mol_%d"%row2
                                if nameCol is not None:
                                    name2 = nameCol.getValueAsString(row2)

                    # print this.selectedColumn.getName(),row1, row2
                                # print strObj1, strObj2
                                if strObj2 is not None:
                                    try:
                                        act2 = float(strObj2)
                                        tanimoto = oegraphsim.oegraphsim.OETanimoto(fp1,fp2)
                                        if tanimoto == 1.0:
                                            tanimoto = 0.9999
                                        sal_score = math.fabs(act1-act2)/(1-tanimoto)
                                        topHits.append((smilesList[row1],name1, act1, smilesList[row2], name2, act2, tanimoto, sal_score))
                                    except:
                                        traceback.print_exc()
                                        continue
                        except:
                            traceback.print_exc()
                            continue
                print "finished"
                hits = sorted(topHits, key=lambda x:x[7], reverse=True)
                n = len(hits)
                if n > this.numHits:
                    n = this.numHits
                return hits[0:n]

            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setNote(pair[0])
                self.progressMonitor.setProgress(pair[1])

            def done(self):
                try:
                    hitsFound = self.get()  #raise exception if abnormal completion
                    vtable2 = VortexTableModel()
                    smiCol1 = vtable2.addColumn("Smiles A",Column.COLUMN_STRING)
                    nameCol1 = vtable2.addColumn("Name A",Column.COLUMN_STRING)
                    actCol1 = vtable2.addColumn("Act A",Column.COLUMN_DOUBLE)
                    smiCol2 = vtable2.addColumn("Smiles B",Column.COLUMN_STRING)
                    nameCol2 = vtable2.addColumn("Name B",Column.COLUMN_STRING)
                    actCol2 = vtable2.addColumn("Act B",Column.COLUMN_DOUBLE)
                    sal_scoreCol = vtable2.addColumn("Activity Cliff Score",Column.COLUMN_DOUBLE)
                    tanimotoCol = vtable2.addColumn("Tanimoto",Column.COLUMN_DOUBLE)
                    for id,obj in enumerate(hitsFound):
                        smiles1,name1,act1,smiles2,name2,act2,tanimoto,sal_score = obj
                        smiCol1.setValueFromString(id,smiles1)
                        nameCol1.setValueFromString(id,name1)
                        actCol1.setValueFromString(id,"%f"%round(act1,2))
                        smiCol2.setValueFromString(id,smiles2)
                        nameCol2.setValueFromString(id,name2)
                        actCol2.setValueFromString(id,"%f"%round(act2,2))
                        sal_scoreCol.setValueFromString(id,"%f"%round(sal_score,2))
                        tanimotoCol.setValueFromString(id,"%f"%round(tanimoto,2))
                    table2 = VortexJTable2()
                    table2.setModel(vtable2)
                    ws = vortex.createWorkspace(vtable2,"Activity Cliff")
                    ws.addComponent(JScrollPane(table2),"Activity Cliff")
                    ws.updateLayout(VortexWorkspace.LAYOUT_GRID)

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
