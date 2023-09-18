#Imports of libraries
import java.lang.String as String
#import javax.xml
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.util.Util as Util
import sys,os

sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

from java.util.concurrent import ExecutionException

from javax.swing.border import TitledBorder
from javax.swing import ProgressMonitor,SwingWorker
from javax.swing import JComboBox
import openeye.oechem.OEGraphMol as OEGraphMol
import openeye.oechem as oechem

class DuplicateTagger:
    def __init__(self):
        return


    def show_interface(self):
        mypanel = swing.JPanel()
        mypanel.setBorder(TitledBorder("Select Structure Column"))


        smi_cols = []
        for column in vtable.getColumnsArray():
            colName = String(column.getName())
            if colName.contains("Smiles") or colName.contains("SMILES") or colName.contains("smiles") or colName.equalsIgnoreCase("Smiles") or colName.equals("Structure"):
                smi_cols.append(colName)
        self.comboBox = JComboBox(smi_cols)

        mypanel.add(self.comboBox)

        ret = vortex.showInDialog(mypanel, "Tag duplicate structuers")

        if ret == vortex.OK:
            self.colName = self.comboBox.getSelectedItem()
            return len(self.colName)>0
        else:
            return False


    def calculateAll(self):
        colName = self.colName
        isStructure = False
        if colName == "Structure":
            isStructure = True

        class MyTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setMillisToDecideToPopup(0)
                self.progressMonitor.setMillisToPopup(0)

            def doInBackground(self):
                useSelected = False
                nrow = vtable.getRealRowCount()
                rows = range(0,nrow)
                dup_column = vtable.findColumnWithName("DuplicateOrNot",1)
                duplicateDB = {}

                for row in rows:
                    if self.progressMonitor.isCanceled():
                        break

                    self.super__publish([(100*row/len(rows),"process compound %d"%(row))])
                    smi_column = vtable.findColumnWithName(colName,0)

                    mol = OEGraphMol()
                    absSmiles = None

                    if isStructure:
                        molfile = vtable.getMolFileManager().getMolFileAtRow(row)
                        oemolistream = oechem.oemolistream()
                        oemolistream.SetFormat(oechem.OEFormat.MDL)
                        oemolistream.openstring(molfile)
                        oechem.oechem.OEReadMolecule(oemolistream,mol)
                        absSmiles = oechem.oechem.OEMolToSmiles(mol)
                    else:
                        smiles = smi_column.getValueAsString(row)
                        oemolistream = oechem.oemolistream()
                        oemolistream.SetFormat(oechem.OEFormat.SMI)
                        oemolistream.openstring(smiles)
                        oechem.oechem.OEReadMolecule(oemolistream,mol)
                        absSmiles = oechem.oechem.OEMolToSmiles(mol)

                    if not duplicateDB.has_key(absSmiles):
                        dup_column.setValueFromString(row,"NO")
                        duplicateDB[absSmiles] = row
                    else:
                        dup_column.setValueFromString(row,"YES")

            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setProgress(pair[0])
                self.progressMonitor.setNote(pair[1])


            def done(self):
                self.progressMonitor.close()
                try:
                    self.get()  #raise exception if abnormal completion
                    vtable.fireTableStructureChanged()
                except ExecutionException, e:
                    e.printStackTrace()
                    vortex.alert(e.getMessage())

        task = MyTask()
        task.execute()

if __name__ == "__main__":
    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        tagger = DuplicateTagger()
        if tagger.show_interface():
            tagger.calculateAll()
