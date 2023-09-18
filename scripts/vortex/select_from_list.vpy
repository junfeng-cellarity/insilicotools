__doc__ = "Plot ClogP vs pKi with lines at LIPE levels."
#Imports of libraries
import os,sys
import java.lang as lang
import java.lang.String as String
import java.net.URL as URL
import java.net.URLConnection as URLConnection
import java.io.InputStreamReader as InputStreamReader
import java.io.BufferedReader as BufferedReader
import java.io.DataOutputStream as DataOutputStream
import java.util.Properties as Properties
import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.util.Util as Util
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.table.VortexJTable2 as VortexJTable2
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.IntArray as IntArray
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.plots.VortexScatterPlot as VortexScatterPlot
import com.dotmatics.vortex.graphicalobjects as grob
import com.dotmatics.vortex.plots.VortexScatterPlot as vp
import traceback
import re

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
from javax.swing import JButton, JFrame, JTextField,JComboBox,JTextPane
from java.awt import GridBagLayout as awtGridBagLayout, BorderLayout, GridLayout, FlowLayout
from java.awt import GridBagConstraints as awtGridBagConstraints
from java.awt import Insets as awtInsets
import java.awt.Color as Color
import java.awt.Font as Font
import java.awt.Dimension as Dimension
import base64
from java.lang import Class
import java.sql.DriverManager as DriverManager


SINGLE_BATCH = 1
SMALL_BATCH = 4
MEDIUM_BATCH = 11
LARGE_BATCH = 51

class BioNumberConversion:
    def __init__(self):
        self.columns = []
        for column in vtable.getColumnsArray():
            # if not (column.isStructureColumn()):
            self.columns.append(column.getName())
        self.columnCB = JComboBox(self.columns)
        self.textpane = JTextPane()

    def show_interface(self):
        if len(self.columns)==0:
            vortex.alert("No suitable column available to convert.")
            return False

        mypanel = JPanel(BorderLayout())


        columnPanel = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        columnLabel = JLabel("Use Column:")
        columnPanel.add(columnLabel)
        columnPanel.add(self.columnCB)

        mypanel.add(JScrollPane(self.textpane),BorderLayout.CENTER)
        mypanel.add(columnPanel,BorderLayout.SOUTH)
        mypanel.setPreferredSize(Dimension(800,600))

        ret = vortex.showInDialog(mypanel, "Select Rows by Column Values:")
        if ret == vortex.OK:
            return True
        else:
            return False

    def calculateAll(self):
        
        columnName = self.columnCB.getSelectedItem()
        buffer = str(self.textpane.getText()).split("\n")

        class ConvertTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setNote("Running ...")
                self.progressMonitor.setProgress(1)
                self.selectedColumn = vtable.findColumnWithName(columnName,0)

            def doInBackground(self):
                self.super__publish([(1,"Running ...")])
                nrow = vtable.getRealRowCount()
                selected_list = []
                for v in buffer:
                    selected_list.append(v.strip())

                for row in range(0,nrow):
                    batch_id = self.selectedColumn.getValueAsString(row).strip()
                    if batch_id in selected_list:
                        vtable.setRowSelected(row,True)
                        self.super__publish([(row*100/nrow,"Running ...")])

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
                    vortex.alert(e.getMessage())

        task = ConvertTask()
        task.execute()
        
if __name__ == "__main__":
    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        lep = BioNumberConversion()
        if lep.show_interface():
            lep.calculateAll()
