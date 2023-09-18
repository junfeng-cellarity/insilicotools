import com.dotmatics.vortex
import com.dotmatics.vortex.util.Util as Util
import com.dotmatics.vortex.table.Column as Column
import sys,os
import base64

sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()


from openeye.oechem import *
from openeye.oedepict import *
from javax.swing.border import TitledBorder
from javax.swing import JLabel
from javax.swing import JComboBox
from javax.swing import ProgressMonitor
from javax.swing import BoxLayout
from java.util.concurrent import ExecutionException
import javax.swing.SwingWorker as SwingWorker
from java.io import File

import com.scitegic.proxy.ComponentDatabase
import com.scitegic.proxy.Job
import com.scitegic.proxy.JobStatus
import com.scitegic.proxy.PipelinePilotServer as PipelinePilotServer
import com.scitegic.proxy.PipelinePilotServerConfig
import com.scitegic.proxy.XmldbItem

package_directory = "Protocols/Discngine/TIBCO Spotfire Enabled Protocols/Biogen/Data Functions/"
projects = [
    "BTK",
    "GSK3ab",
    "cnsIRAK4",
    "OGA"
]
protocols_dict = {
    "BTK":"DF SAR Analysis BTK",
    "GSK3ab":"DF SAR Analysis GSK3ab",
    "cnsIRAK4":"DF SAR Analysis cnsIRAK4",
    "OGA":"DF SAR Analysis OGA"
}

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
        oedepict.OEAddDepictionHydrogens(mol)
    else:
        ifs = oemolistream()
        ifs.SetFormat(OEFormat.SDF)
        ifs.openstring(molString)
        oechem.OEReadMolecule(ifs, mol)
        oechem.OESuppressHydrogens(mol)

    if mol.GetDimension()!=2:
        oedepict.OEDepictCoordinates(mol)
        oechem.OEMDLPerceiveBondStereo(mol)

    mol.SetTitle(name)
    if mol.NumAtoms() == 0:
        print >> sys.stderr, "Unable to parse molecule: %s"%molString
        #raise Exception("Unable to parse molecule: %s"%molString)
    return (mol)


class DataHandler:
    def __init__(self):
        self.isCancelled = False
        self.protocol = None
        return


    def show_interface(self):
        mypanel = swing.JPanel()
        layout = BoxLayout(mypanel,BoxLayout.Y_AXIS)
        mypanel.setLayout(layout)

        mypanel.setBorder(TitledBorder("Projects"))
        p1 = swing.JPanel()
        mypanel.add(p1)

        projectsCB = JComboBox(projects)

        p1.add(JLabel("Select Project to Run:"))
        p1.add(projectsCB)

        ret = vortex.showInDialog(mypanel, "Run Pipeline Pilot Protocol:")
        if ret == vortex.OK:
            self.isCancelled = False
            project = projectsCB.getSelectedItem()
            filename = protocols_dict[project]
            self.protocol = "%s/%s/%s"%(package_directory,project,filename)
            return len(self.protocol)>0
        else:
            self.isCancelled = True
            return False


    def calculateAll(self):
        this = self
        class MyTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setMillisToDecideToPopup(0)
                self.progressMonitor.setMillisToPopup(0)

            def doInBackground(self):
                self.super__publish([(1,"Running ...")])
                pp = PipelinePilotServer("pipelpp01:9944","pilotweb","pilotweb")
                job = pp.createJob(this.protocol)
                result = job.runAndPoll()
                status = job.getStatus()
                if status.Finished.equals(status):
                    results = result.getResultFiles()
                    if len(results) > 0:
                        for r in results:
                            tmpSdf = com.dotmatics.vortex.Vortex.getTempFile('sdf')
                            pp.getRemoteFileManager().downloadFile(r,File(tmpSdf))
                t = com.dotmatics.vortex.table.VortexTableModel.readSDF(tmpSdf)
                return os.path.basename(r),t

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

