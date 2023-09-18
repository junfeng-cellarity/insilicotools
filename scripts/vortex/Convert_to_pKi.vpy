#Imports of libraries
import java.lang as lang
import java.lang.String as String
import java.net.URL as URL
import java.net.URLConnection as URLConnection
import java.io.InputStreamReader as InputStreamReader
import java.io.BufferedReader as BufferedReader
import java.io.DataOutputStream as DataOutputStream
import java.util.Properties as Properties
#import javax.xml
import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.util.Util as Util
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.table.VortexJTable2 as VortexJTable2
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.IntArray as IntArray
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.plots.VortexScatterPlot as VortexScatterPlot

__doc__ = 'Convert Ki or IC50 to pKi or pIC50'
__author__ = 'Jun Feng'
__version__ = '1.0'

import sys
import time
import java.net.URL
import math
import jarray
import traceback

from java.lang import Runnable, Double, Integer, Boolean
from java.util import Random, HashMap
from java.util.concurrent import ExecutionException

from javax.swing import BorderFactory, ProgressMonitor, SwingWorker, JScrollPane, JComboBox, JPanel, JLabel
from java.text import DecimalFormat
from javax.swing import JButton, JFrame, JTextField
from java.awt import GridBagLayout as awtGridBagLayout, BorderLayout, GridLayout, FlowLayout
from java.awt import GridBagConstraints as awtGridBagConstraints
from java.awt import Insets as awtInsets

class LogAct:
    def __init__(self):
        self.concDict = {"nM":1e-9,"uM":1e-6,"mM":1e-3,"None":1}
        self.modeCB = JComboBox(["-Log","Log"])
        self.concCB = JComboBox(["nM","uM","mM","None"])
        self.concCB.setSelectedItem("uM")
        self.columns = []
        for column in vtable.getColumnsArray():
            if column.getForceNumeric() or column.isNumeric():
                self.columns.append(column.getName())
        self.columnCB = JComboBox(self.columns)

    def show_interface(self):
        if len(self.columns)==0:
            vortex.alert("No numeric column available to convert.")
            return False

        mypanel = JPanel(GridLayout(3,1))

        concPanel = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        concLabel = JLabel("Concentration:")
        concPanel.add(concLabel)
        concPanel.add(self.concCB)

        modePanel = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        modeLabel = JLabel("Function:")
        modePanel.add(modeLabel)
        modePanel.add(self.modeCB)

        columnPanel = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        columnLabel = JLabel("Use Column:")
        columnPanel.add(columnLabel)
        columnPanel.add(self.columnCB)

        mypanel.add(columnPanel)
        mypanel.add(modePanel)
        mypanel.add(concPanel)

        ret = vortex.showInDialog(mypanel, "Convert to pKi")
        if ret == vortex.OK:
            return True
        else:
            return False


    def calculateAll(self):

        columnName = self.columnCB.getSelectedItem()
        functionName = String(self.modeCB.getSelectedItem())
        conc = self.concCB.getSelectedItem()
        convertRatio = self.concDict[conc]

        class ConvertTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setNote("Running ...")
                self.progressMonitor.setProgress(1)
                self.selectedColumn = vtable.findColumnWithName(columnName,False)
                self.newColumn = vtable.findColumnWithName("%s_%s(%s)"%(functionName,columnName,conc),True,Column.COLUMN_DOUBLE)

            def doInBackground(self):
                nrow = vtable.getRealRowCount()
                for row in range(0, nrow):
                    if self.progressMonitor.isCanceled():
                        break
                    self.super__publish([(100*row/nrow,"process compound %d"%(row))])
                    try:
#                    if 1:
                        columnValue = String(self.selectedColumn.getValueAsString(row))
                        columnValue2 = String(columnValue.replaceAll(">", ""))
                        originValue = float(Double.parseDouble(columnValue2.replaceAll("<","")))
                        originValue = math.log10(originValue*convertRatio)
                        if not functionName.equals("Log"):
                            originValue = -originValue
                        self.newColumn.setDouble(row,round(originValue,2))
                    except:
                        print >>sys.stderr, self.selectedColumn.getValueAsString(row)
                        self.newColumn.setValueFromString(row,None)

            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setProgress(pair[0])
                self.progressMonitor.setNote(pair[1])


            def done(self):
                self.progressMonitor.close()
                try:
                    self.get()  #raise exception if abnormal completion
                    vtable.fireTableStructureChanged()
                    vortex.alert("A new column called %s has been added to your Vortex table."%(self.newColumn.getName()))
                except ExecutionException, e:
                    vortex.alert(e.getMessage())

        task = ConvertTask()
        task.execute()

if __name__ == "__main__":
    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        logAct = LogAct()
        if logAct.show_interface():
            logAct.calculateAll()