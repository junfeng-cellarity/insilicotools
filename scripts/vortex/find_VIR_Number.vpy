__doc__ = ""
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
import urllib, subprocess, xmlrpclib


dbURL = "jdbc:postgresql://10.74.2.128/skyhawk_virtual"
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

PORT = 9527
MODELING_SERVER = "http://10.74.2.128:%d/" % PORT
class ModelingServer:
    def __init__(self):
        self.server =  xmlrpclib.Server(MODELING_SERVER)

    def get_canonicalized_structure(self,mol_file):
        return self.server.canonicalizer(mol_file)

class RGXNumberSearch:
    def __init__(self):
        self.columns = []
        for column in vtable.getColumnsArray():
            if column.isStructureColumn():
                self.columns.append(column.getName())
        self.columnCB = JComboBox(self.columns)
        self.server = ModelingServer()

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

        ret = vortex.showInDialog(mypanel, "Pick Structure Column:")
        if ret == vortex.OK:
            return True
        else:
            return False

    def calculateAll(self):
        
        columnName = self.columnCB.getSelectedItem()
        this = self
        class ConvertTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                Class.forName("oracle.jdbc.OracleDriver")
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setNote("Running ...")
                self.progressMonitor.setProgress(1)
                # self.selectedColumn = vtable.findColumnWithName(columnName,False)

            def getRGXNumberFromMol(self,molfile):
                rs = None
                statement = None
                connection = None
                try:
                    connection = DriverManager.getConnection(dbURL, db_user, db_password)
                    statement = connection.prepareStatement("select corporate_id from chemicals a inner join rdk.mols b on a.chemical_id=b.chemical_id where m@=mol_from_ctab(?::cstring)")
                    statement.setString(1,molfile)
                    rs = statement.executeQuery()
                    if rs.next():
                        rgx_number = rs.getString("corporate_id")
                        return rgx_number
                    else:
                        return None
                finally:
                    if rs is not None:
                        rs.close()
                    if statement is not None:
                        statement.close()
                    if connection is not None:
                        connection.close()


            # def getRGXNumberFromSmiles(self,smiles):
            #     rs = None
            #     statement = None
            #     connection = None
            #     try:
            #         connection = DriverManager.getConnection(dbURL, db_user, db_password)
            #         statement = connection.prepareStatement("select corporate_id from chemicals a inner join rdk.mols b on a.chemical_id=b.chemical_id where m@=mol_from_smiles(?::cstring)")
            #         statement.setString(1,smiles)
            #         rs = statement.executeQuery()
            #         if rs.next():
            #             rgx_number = rs.getString("corporate_id")
            #             return rgx_number
            #         else:
            #             return None
            #     finally:
            #         if rs is not None:
            #             rs.close()
            #         if statement is not None:
            #             statement.close()
            #         if connection is not None:
            #             connection.close()


            def doInBackground(self):
                self.super__publish([(1,"Running ...")])
                nrow = vtable.getRealRowCount()
                col = None
                new_col = vtable.findColumnWithName("VIR_Number",True,Column.COLUMN_STRING)
                if columnName!="Structure":
                    col = vtable.findColumnWithName(columnName,False)

                for row in range(0, nrow):
                    self.super__publish([(100*(row)/nrow,"process compound %d"%(row))])
                    try:
                        rgx_number = None
                        if columnName == "Structure":
                            mfm = vtable.getMolFileManager()
                            molfile = mfm.getMolFileAtRow(row)
                            oemol = OEGraphMol()
                            oechem.OEReadMolFromString(oemol,OEFormat.MDL,False,molfile)
                            oechem.OEDeleteEverythingExceptTheFirstLargestComponent(oemol)
                            molfile = oechem.OEWriteMolToString(OEFormat.MDL,False, oemol)
                            molfile = this.server.get_canonicalized_structure(molfile)
                            rgx_number = self.getRGXNumberFromMol(molfile)
                        else:
                            oemol = OEGraphMol()
                            smiles = col.getValueAsString(row)
                            oechem.OEParseSmiles(oemol,smiles)
                            oechem.OEDeleteEverythingExceptTheFirstLargestComponent(oemol)
                            oechem.OESuppressHydrogens(oemol)
                            molfile = oechem.OEWriteMolToString(OEFormat.MDL,False, oemol)
                            molfile = this.server.get_canonicalized_structure(molfile)
                            rgx_number = self.getRGXNumberFromMol(molfile)
                        if rgx_number is not None:
                            new_col.setValueFromString(row,rgx_number)
                    except:
                        traceback.print_exc()

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
        lep = RGXNumberSearch()
        if lep.show_interface():
            lep.calculateAll()
