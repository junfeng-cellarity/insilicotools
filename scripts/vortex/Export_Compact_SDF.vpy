__doc__ = "Extract a subset of molecules, remove empty columns and export to a new workspace"
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

sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()
###

from chemaxon.formats import *
from chemaxon.struc import *

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
from javax.swing import JButton, JFrame, JTextField,JComboBox
from java.awt import GridBagLayout as awtGridBagLayout, BorderLayout, GridLayout, FlowLayout
from java.awt import GridBagConstraints as awtGridBagConstraints
from java.awt import Insets as awtInsets
import java.awt.Color as Color
import java.awt.Font as Font
import base64
from java.lang import Class
import java.sql.DriverManager as DriverManager

import openeye.oechem.OEFormat as OEFormat
import openeye.oechem.OEGraphMol as OEGraphMol
import openeye.oechem.oemolistream as oemolistream
import openeye.oechem.oemolostream as oemolostream
from openeye.oechem import *
from openeye.oedepict import *


def convertToOEMol(molString, name):
    mol = OEGraphMol()
    ifs = oemolistream()
    ifs.SetFormat(OEFormat.SDF)
    ifs.openstring(molString)
    oechem.OEReadMolecule(ifs, mol)
    mol.SetTitle(name)
    if mol.NumAtoms() == 0:
        print >> sys.stderr, "Unable to parse molecule: %s"%molString
        #raise Exception("Unable to parse molecule: %s"%molString)
    return (mol)


class CompactSDFConversion:
    def __init__(self):
        return

    def calculateAll(self):
        
        class ConvertTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setNote("Running ...")
                self.progressMonitor.setProgress(1)

            def doInBackground(self):
                self.super__publish([(1,"Running ...")])
                selected_rows = vtable.getSelectedRows()
                selected_columns = []
                for column in vtable.getColumnsArray():
                    if column.getName().strip()=="Structure":
                        continue
                    isEmpty = True
                    for row_id in range(0, vtable.getSelectedRowCount()):
                        row = selected_rows.get(row_id)
                        col_value = column.getValueAsString(row).strip()
                        if col_value is not None and len(col_value) > 0 and col_value != "N/A":
                            isEmpty = False
                            break
                    if not isEmpty:
                        selected_columns.append(column)
                newName = vtable.name + ".sdf" if vtable.name.find(".") == -1 else vtable.name[:vtable.name.find(".")] + ".sdf"
                tmpSdf = com.dotmatics.vortex.Vortex.getTempFile('sdf')

                ofs = oemolostream(tmpSdf)
                for row_id in range(0, vtable.getSelectedRowCount()):
                    row = selected_rows.get(row_id)
                    try:
                        mfm = vtable.getMolFileManager()
                        molfile = mfm.getMolFileAtRow(row)
                        oemol = OEGraphMol()
                        oechem.OEReadMolFromString(oemol,OEFormat.MDL,False,molfile)
                        for col in selected_columns:
                            oechem.OESetSDData(oemol,col.getName().encode("utf-8"),col.getValueAsString(row).encode("utf-8"))
                        oechem.OEWriteMolecule(ofs,oemol)

                    except:
                        pass
                    #     print >>sys.stderr, row
                    #     traceback.print_exc()
                    #     oemol = OEGraphMol()
                    #     for c in range(0,vtable.getColumnCount()):
                    #         name = vtable.getColumnName(c)
                    #         col = vtable.getColumn(c)
                    #         oechem.OESetSDData(oemol,name,col.getValueAsString(row))
                    #     oechem.OEWriteMolecule(ofs,oemol)
                ofs.close()
                t = com.dotmatics.vortex.table.VortexTableModel.readSDF(tmpSdf)
                vortex.addTable(newName, t)


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
        lep = CompactSDFConversion()
        lep.calculateAll()
