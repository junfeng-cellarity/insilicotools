import com.dotmatics.vortex
import com.dotmatics.vortex.util.Util as Util
import com.dotmatics.vortex.table.Column as Column
import sys,os
import urllib2
import time

# sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
# from initOEChem import VortexInit
#
# vortexInit = VortexInit(vortex)
# vortexInit.initOpeneye()
# vortexInit.initJars()


# from openeye.oechem import *
# from openeye.oedepict import *
from javax.swing.border import TitledBorder
from javax.swing import JLabel
from javax.swing import JComboBox
from javax.swing import ProgressMonitor
from javax.swing import BoxLayout
from java.util.concurrent import ExecutionException
import javax.swing.SwingWorker as SwingWorker
from java.io import File

import javax.swing.SwingUtilities as SwingUtilities

prj_url = "http://10.74.2.128:8080/data/%s_project.sdf"
tplt_url = "http://10.74.2.128:8080/data/%s_template.vortexgz"
projects = [
    "htt",
    "sma",
    "cMYB",
#    "ar",
#    "wrn",
#    "covid",
    "mdck",
    "hERG"
]

class DataHandler:
    def __init__(self):
        self.isCancelled = False
        self.project_url = None
        self.project_name = None
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

        ret = vortex.showInDialog(mypanel, "Load Project Data")
        if ret == vortex.OK:
            self.isCancelled = False
            self.project_name = projectsCB.getSelectedItem()
            self.project_url = prj_url % (self.project_name)
            self.template_url = tplt_url%(self.project_name)
            return (self.project_url is not None)
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
                response = urllib2.urlopen(this.project_url)
                if this.project_url.endswith("sdf"):
                    sdf = response.read()
                    sdf = sdf.replace('\r\n', '\n')
                    tmpSdf = com.dotmatics.vortex.Vortex.getTempFile('sdf')
                    binary_io = open(tmpSdf, "w")
                    binary_io.write(sdf)
                    binary_io.close()
                    vortexgzFile = None
                    try:
                        response = urllib2.urlopen(this.template_url)
                        vortexgz = response.read()
                        vortexgzFile = com.dotmatics.vortex.Vortex.getTempFile("vortexgz")
                        b = open(vortexgzFile, "wb")
                        b.write(vortexgz)
                        b.close()
                    except Exception,ex:
                        print(ex)
                    return this.project_name,tmpSdf,vortexgzFile
                else:
                    csv = response.read()
                    tmp_csv = com.dotmatics.vortex.Vortex.getTempFile('csv')
                    binary_io = open(tmp_csv, "w")
                    binary_io.write(csv)
                    binary_io.close()
                    return this.project_name,tmp_csv,None


            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setProgress(pair[0])
                self.progressMonitor.setNote(pair[1])


            def done(self):
                try:
                    newName,sdf,template = self.get()
                    if template is not None:
                        table = vortex.readSDF(sdf)
                        ws = vortex.createWorkspace(table, newName)
                        ws.applyTemplate(template)
                        ws.updateLayout(VortexWorkspace.LAYOUT_GRID)
                    else:
                        vortex.loadAnyFile(sdf)
                finally:
                    self.progressMonitor.close()


        task = MyTask()
        task.execute()

if __name__ == "__main__":

    dataHandler = DataHandler()

    while not dataHandler.isCancelled:
        if dataHandler.show_interface():
            break

    if not dataHandler.isCancelled:
        dataHandler.calculateAll()

