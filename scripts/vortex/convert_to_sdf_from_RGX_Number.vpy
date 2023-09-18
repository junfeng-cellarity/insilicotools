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


dbURL = "jdbc:postgresql://10.74.2.128/skyhawk"
db_user = "medchem"
db_password = "medchem"

SINGLE_BATCH = 1
SMALL_BATCH = 4
MEDIUM_BATCH = 11
LARGE_BATCH = 51

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


class RGXNumberConversion:
    def __init__(self):
        self.columns = []
        for column in vtable.getColumnsArray():
            if not (column.getForceNumeric() or column.isNumeric() or column.isStructureColumn()):
                self.columns.append(column.getName())
        self.columnCB = JComboBox(self.columns)

    def show_interface(self):
        if len(self.columns)==0:
            vortex.alert("No suitable column available to convert.")
            return False

        mypanel = JPanel()

        columnPanel = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        columnLabel = JLabel("Use Column:")
        columnPanel.add(columnLabel)
        columnPanel.add(self.columnCB)

        mypanel.add(columnPanel)

        ret = vortex.showInDialog(mypanel, "Pick Column with BioNumber:")
        if ret == vortex.OK:
            return True
        else:
            return False

    def calculateAll(self):
        
        columnName = self.columnCB.getSelectedItem()

        class ConvertTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                Class.forName("oracle.jdbc.OracleDriver")
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setNote("Running ...")
                self.progressMonitor.setProgress(1)
                self.selectedColumn = vtable.findColumnWithName(columnName,False)

            def getMolListFromRGXNumberListBatch(self, rgx_number_list):
                if rgx_number_list is None or len(rgx_number_list) == 0:
                    return []
                statements = []
                resultsets = []
                molFileDict = {}
                connection = None
                try:
                    connection = DriverManager.getConnection(dbURL, db_user, db_password)
                    totalNumberOfValuesLeftToBatch = len(rgx_number_list)
                    total = len(rgx_number_list)
                    while totalNumberOfValuesLeftToBatch>0:
                        batch_size = SINGLE_BATCH
                        if totalNumberOfValuesLeftToBatch>= LARGE_BATCH:
                            batch_size = LARGE_BATCH
                        else:
                            if totalNumberOfValuesLeftToBatch >= MEDIUM_BATCH:
                                batch_size = MEDIUM_BATCH
                            else:
                                if totalNumberOfValuesLeftToBatch >= SMALL_BATCH:
                                    batch_size = SMALL_BATCH
                        totalNumberOfValuesLeftToBatch -= batch_size
                        self.super__publish([(100*(total-totalNumberOfValuesLeftToBatch)/total,"process compound %d"%(total-totalNumberOfValuesLeftToBatch))])
                        tmplist1 =  []
                        for i in range(0,batch_size):
                            tmplist1.append("?")
                        clause = ",".join(tmplist1)
                        statement = connection.prepareStatement("select orig_molfile, corporate_id from chemicals where corporate_id in (%s)"%clause)
                        for i in range(1,batch_size+1):
                            rgx_number = rgx_number_list[i-1]
                            statement.setString(i,rgx_number)

                        for i in range(batch_size-1,-1,-1):
                            del(rgx_number_list[i])
                        statements.append(statement)
                        rs = statement.executeQuery()
                        resultsets.append(rs)
                        while rs.next():
                            molFile = rs.getString("orig_molfile")
                            rgx_number = rs.getString("corporate_id")
                            molFileDict[rgx_number] = molFile
                        try:
                            rs.close()
                            statement.close()
                        except:
                            pass
                    return molFileDict


                finally:
                    for rs in resultsets:
                        if rs is not None:
                            try:
                                rs.close()
                            except:
                                traceback.print_exc()
                    for ps in statements:
                        if ps is not None:
                            try:
                                ps.close()
                            except:
                                traceback.print_exc()
                    if connection is not None:
                        try:
                            connection.close()
                        except:
                            traceback.print_exc()

            def doInBackground(self):
                self.super__publish([(1,"Running ...")])
                nrow = vtable.getRealRowCount()
                rgx_number_list = []
                for row in range(0, nrow):
                    try:
                        rgx_number = str(self.selectedColumn.getValueAsString(row))
                        rgx_number = re.sub(r"RGX-([0-9][0-9][0-9][0-9][0-9][0-9][0-9])-([1-9])",r"RGX-\1",rgx_number)
                        rgx_number = re.sub(r"RGX-([0-9][0-9][0-9][0-9][0-9][0-9][0-9])-(0[1-9])",r"RGX-\1",rgx_number)
                        rgx_number_list.append(rgx_number)
                    except:
                        print >>sys.stderr, self.selectedColumn.getValueAsString(row)
                        traceback.print_exc()
                molDict = self.getMolListFromRGXNumberListBatch(rgx_number_list)

                newName = vtable.name + ".sdf" if vtable.name.find(".") == -1 else vtable.name[:vtable.name.find(".")] + ".sdf"
                tmpSdf = com.dotmatics.vortex.Vortex.getTempFile('sdf')

                ofs = oemolostream(tmpSdf)
                for row in range(0,nrow):
                    batch_id = self.selectedColumn.getValueAsString(row)
                    try:
                        rgx_number = str(batch_id)
                        rgx_number = re.sub(r"RGX-([0-9][0-9][0-9][0-9][0-9][0-9][0-9])-([1-9])",r"RGX-\1",rgx_number)
                        rgx_number = re.sub(r"RGX-([0-9][0-9][0-9][0-9][0-9][0-9][0-9])-(0[1-9])",r"RGX-\1",rgx_number)
                        molFile = molDict[rgx_number]
                        oemol = convertToOEMol(molFile, batch_id)
                        for c in range(0,vtable.getColumnCount()):
                            name = vtable.getColumnName(c)
                            col = vtable.getColumn(c)
                            oechem.OESetSDData(oemol,name,col.getValueAsString(row))
                        oechem.OEDeleteEverythingExceptTheFirstLargestComponent(oemol)
                        oechem.OEWriteMolecule(ofs,oemol)

                    except:
                        print >>sys.stderr, batch_id
                        traceback.print_exc()
                        oemol = OEGraphMol()
                        for c in range(0,vtable.getColumnCount()):
                            name = vtable.getColumnName(c)
                            col = vtable.getColumn(c)
                            oechem.OESetSDData(oemol,name,col.getValueAsString(row))

                        oechem.OEWriteMolecule(ofs,oemol)
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
        lep = RGXNumberConversion()
        if lep.show_interface():
            lep.calculateAll()
