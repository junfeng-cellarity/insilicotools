import com.dotmatics.vortex
import com.dotmatics.vortex.util.Util as Util
import com.dotmatics.vortex.table.Column as Column
import sys,os

sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()


from openeye.oechem import *
from javax.swing.border import TitledBorder
from javax.swing import JLabel
from javax.swing import JComboBox
from javax.swing import ProgressMonitor
from javax.swing import BoxLayout
from java.util.concurrent import ExecutionException
import javax.swing.SwingWorker as SwingWorker

def getSmilesOrMolColumns(vtable):
    structure_cols = []
    for column in vtable.getColumnsArray():
        colName = column.getName()
        lowerColName = colName.lower()
        if "smiles" in lowerColName:
            structure_cols.append(colName)
        if "molfile" in lowerColName:
            structure_cols.append(colName)
    return structure_cols


def convertToOEMol(molString, name, isSmiles):
    mol = OEGraphMol()
    if isSmiles:
        oechem.OEParseSmiles(mol, molString)
        oechem.OESuppressHydrogens(mol)
        oechem.OEGenerate2DCoordinates(mol)
    else:
        ifs = oemolistream()
        ifs.SetFormat(OEFormat.SDF)
        ifs.openstring(molString)
        oechem.OEReadMolecule(ifs, mol)
        oechem.OESuppressHydrogens(mol)

    if mol.GetDimension()!=2:
        oechem.OEMDLPerceiveBondStereo(mol)
        oechem.OEGenerate2DCoordinates(mol)

    mol.SetTitle(name)
    if mol.NumAtoms() == 0:
        print >> sys.stderr, "Unable to parse molecule: %s"%molString
        #raise Exception("Unable to parse molecule: %s"%molString)
    return (mol)


class DataHandler:
    def __init__(self):
        self.isCancelled = False
        return


    def show_interface(self):
        structure_cols = getSmilesOrMolColumns(vtable)

        if len(structure_cols) > 0:
            mypanel = swing.JPanel()
            layout = BoxLayout(mypanel,BoxLayout.Y_AXIS)
            mypanel.setLayout(layout)

            mypanel.setBorder(TitledBorder("Columns"))
            p1 = swing.JPanel()
            p2 = swing.JPanel()
            mypanel.add(p1)
            mypanel.add(p2)

            tag_cols = []
            for column in vtable.getColumnsArray():
                colName = java.lang.String(column.getName())
                tag_cols.append(colName)
            tagCB = JComboBox(tag_cols)

            p1.add(JLabel("Select Name Column:"))
            p1.add(tagCB)

            self.scb = JComboBox(structure_cols)
            p2.add(JLabel("Select Smiles column:"))
            p2.add(self.scb)
            ret = vortex.showInDialog(mypanel, "Convert to SDF")
            if ret == vortex.OK:
                self.isCancelled = False
                self.tag = tagCB.getSelectedItem()
                return len(self.tag)>0
            else:
                self.isCancelled = True
                return False
        else:
            vortex.alert("No Smiles column available.")
            self.isCancelled = True
            return False


    def calculateAll(self):
        tag = self.tag
        structural_col = self.scb.getSelectedItem()
        lowerColName = structural_col.lower()
        if "smiles" in lowerColName:
            useSmiles = True
        if "molfile" in lowerColName:
            useSmiles = False

        class MyTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setMillisToDecideToPopup(0)
                self.progressMonitor.setMillisToPopup(0)

            def doInBackground(self):
                self.super__publish([(1,"Running ...")])
                (structureCol, nameCol, tblCols) = (vtable.findColumnWithName(structural_col,0), None, [])

                nameCol = vtable.findColumnWithName(tag,0)
                for (name, col) in [(vtable.getColumnName(c), vtable.getColumn(c)) for c in xrange(vtable.getColumnCount()) if vtable.getColumn(c) != structureCol]:
                    print name
                    if nameCol.getName() == name:
                        continue
                    if "name" in name.lower() and nameCol is None:
                        nameCol = col
                    else:
                        tblCols.append((name, col))
                newName = vtable.name + ".sdf" if vtable.name.find(".") == -1 else vtable.name[:vtable.name.find(".")] + ".sdf"
                tmpSdf = com.dotmatics.vortex.Vortex.getTempFile('sdf')

                ofs = oemolostream(tmpSdf)
                n = 0
                for r in xrange(vtable.getRealRowCount()):
                    self.super__publish([(100*n/vtable.getRealRowCount(),"process mol %d"%n)])
                    name = "unk" if nameCol is None else nameCol.getValueAsString(r)
                    if useSmiles:
                        smi = structureCol.getValueAsString(r)
                        oeMol = convertToOEMol(smi, name,True)
                    else:
                        molfile = structureCol.getValueAsString(r)
                        oeMol = convertToOEMol(molfile,name,False)

                    for (name, col) in tblCols:
                        value = col.getValueAsString(r)
                        if value is not None:
                            oechem.OESetSDData(oeMol, name, value)
                    oechem.OEWriteMolecule(ofs, oeMol)
                    n+=1
                ofs.close()
                t = com.dotmatics.vortex.table.VortexTableModel.readSDF(tmpSdf)
                return newName,t

            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                print pair
                self.progressMonitor.setProgress(pair[0])
                self.progressMonitor.setNote(pair[1])


            def done(self):
                self.progressMonitor.close()
                try:
                    newName,t = self.get()
                    vortex.addTable(newName, t)
                except ExecutionException, e:
                    e.printStackTrace()
                    vortex.alert(e.getMessage())

        task = MyTask()
        task.execute()

if __name__ == "__main__":

    dataHandler = DataHandler()

    while not dataHandler.isCancelled:
        if dataHandler.show_interface():
            break

    if not dataHandler.isCancelled:
        dataHandler.calculateAll()

