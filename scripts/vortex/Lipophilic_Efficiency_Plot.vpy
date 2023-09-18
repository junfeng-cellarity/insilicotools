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

sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()
###

from chemaxon.formats import *
from chemaxon.struc import *
from chemaxon.marvin.plugin.CalculatorPlugin import *
from chemaxon.marvin.calculations import *

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

LOGP_OPTIONS = ["Recalculate LogP","Recalculate LogD7.4"]
class LEPlot:
    def __init__(self):
        self.concDict = {"pM":1e-12,"nM":1e-9,"uM":1e-6,"mM":1e-3,"None":1}
        self.modeCB = JComboBox(["-Log","Log"])
        self.concCB = JComboBox(["pM","nM","uM","mM","None"])
        self.concCB.setSelectedItem("uM")
        self.logpColumns = [LOGP_OPTIONS[0],LOGP_OPTIONS[1]]
        self.columns = []
        for column in vtable.getColumnsArray():
            if column.getForceNumeric() or column.isNumeric():
                self.columns.append(column.getName())
                col_name = str(column.getName()).upper()
                if "LOGP" in col_name or "LOGD" in col_name:
                    self.logpColumns.append(column.getName())

        self.logpColumnCB = JComboBox(self.logpColumns)
        self.columnCB = JComboBox(self.columns)

        self.pkaPlugin = pKaPlugin()
        pkaParams = Properties()
        pkaParams.put("type","PKa,acidic,basic")
        self.pkaPlugin.setParameters(pkaParams)
        self.pkaPlugin.setpH(7.4)
        self.pkaPlugin.setBasicpKaLowerLimit(-5)
        self.pkaPlugin.setAcidicpKaUpperLimit(20)

        self.logPlugin = logPPlugin()

        self.logDPlugin = logDPlugin()
        self.logDPlugin.setpH(7.4)
        
        self.polarizabilityPlugin = PolarizabilityPlugin()
        self.polarizabilityPlugin.setpH(7.4)
        self.refractivityPlugin = RefractivityPlugin()
        
        self.availableProperties = ["CLogP"]
        self.selectedProperties = ["CLogP"]
        

    def show_interface(self):
        if len(self.columns)==0:
            vortex.alert("No numeric column available to convert.")
            return False

        mypanel = JPanel(GridLayout(4,1))

        concPanel = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        concLabel = JLabel("Concentration:")
        concPanel.add(concLabel)
        concPanel.add(self.concCB)

        modePanel = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        modeLabel = JLabel("Function:")
        modePanel.add(modeLabel)
        modePanel.add(self.modeCB)

        columnPanel1 = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        columnLabel1 = JLabel("Acitivty Column:")
        columnPanel1.add(columnLabel1)
        columnPanel1.add(self.columnCB)

        columnPanel2 = JPanel(FlowLayout(FlowLayout.CENTER,20,10))
        columnLabel2 = JLabel("LogP Column:")
        columnPanel2.add(columnLabel2)
        columnPanel2.add(self.logpColumnCB)

        mypanel.add(columnPanel1)
        mypanel.add(columnPanel2)
        mypanel.add(modePanel)
        mypanel.add(concPanel)

        ret = vortex.showInDialog(mypanel, "Convert to pKi")
        if ret == vortex.OK:
            return True
        else:
            return False


    def calculate_cLogD74(self,mol):
        plugin = self.logDPlugin
        plugin.setMolecule(mol)
        plugin.run()
        return plugin.getlogD()

    def calculate_cLogP(self,mol):
        plugin = self.logPlugin
        plugin.setMolecule(mol)
        plugin.run()
        return plugin.getlogPTrue()
    
    def calculateAll(self):
        
        columnName = self.columnCB.getSelectedItem()
        logpColumnName = self.logpColumnCB.getSelectedItem()
        functionName = String(self.modeCB.getSelectedItem())
        conc = self.concCB.getSelectedItem()
        convertRatio = self.concDict[conc]
        
        structure_cols = []
        for column in vtable.getColumnsArray():
            colName = String(column.getName())
            if colName.contains("SMILES") or colName.contains("smiles") or colName.equalsIgnoreCase("Structure"):
                structure_cols.append(colName)
        
        structureCol = structure_cols[0]
        jchem = self

        class ConvertTask(SwingWorker):
            "Class implementing long running task as a SwingWorker thread"

            #######################################################################
            def __init__(self):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setNote("Running ...")
                self.progressMonitor.setProgress(1)
                self.selectedColumn = vtable.findColumnWithName(columnName,False)
                self.newColumn = vtable.findColumnWithName("%s_%s(%s)"%(functionName,columnName,conc),True,Column.COLUMN_DOUBLE)
                self.lipeCol = vtable.findColumnWithName("LiPE",True,Column.COLUMN_DOUBLE)
                

            def doInBackground(self):
                nrow = vtable.getRealRowCount()
                first = True
                logp_column = None
                
                toCal = False
                if logpColumnName == LOGP_OPTIONS[0] or logpColumnName == LOGP_OPTIONS[1]:
                    toCal = True

                pkiList = []
                
                for row in range(0, nrow):
                    if self.progressMonitor.isCanceled():
                        break
                    self.super__publish([(100*row/nrow,"process compound %d"%(row))])
                    try:
                        columnValue = String(self.selectedColumn.getValueAsString(row))
                        columnValue2 = String(columnValue.replaceAll(">", ""))
                        originValue = float(Double.parseDouble(columnValue2.replaceAll("<","")))
                        originValue = math.log10(originValue*convertRatio)
                        if not functionName.equals("Log"):
                            originValue = -originValue
                        pkiList.append(round(originValue,2))
                        self.newColumn.setDouble(row,round(originValue,2))
                        if toCal:
                            if String(structureCol).equalsIgnoreCase("Structure"):
                                molfile = vtable.getMolFileManager().getMolFileAtRow(row)
                            else:
                                structure_column = vtable.findColumnWithName(structureCol,False)
                                molfile = structure_column.getValueAsString(row)
    
                            mol = MolImporter.importMol(molfile)
                            if logpColumnName == LOGP_OPTIONS[0]:
                                if first:
                                    logp_column = vtable.findColumnWithName("CHEMAXON_LOGP",True,Column.COLUMN_DOUBLE)
                                    logp_column.setDecimalPlaces(2)
                                clp = jchem.calculate_cLogP(mol)
                                logp_column.setDouble(row,round(clp,2))
                                self.lipeCol.setDouble(row, round(originValue-clp,2))
                                if first:
                                    first = False
                            else:
                                if first:
                                    logp_column = vtable.findColumnWithName("CHEMAXON_LOGD7.4",True,Column.COLUMN_DOUBLE)
                                    logp_column.setDecimalPlaces(2)
                                clp = jchem.calculate_cLogD74(mol)
                                logp_column.setDouble(row,round(clp,2))
                                self.lipeCol.setDouble(row, round(originValue-clp,2))
                                if first:
                                    first = False
                        else:
                            logp_column = vtable.findColumnWithName(logpColumnName,False)
                            clp = float(logp_column.getValueAsString(row))
                            self.lipeCol.setDouble(row, round(originValue-clp,2))

                    except:
                        print >>sys.stderr, self.selectedColumn.getValueAsString(row)
                        self.newColumn.setValueFromString(row,None)
                clogPList = []
                for row in xrange(0,nrow):
                    clogPList.append(logp_column.getValue(row))

                logp_column_name = logp_column.getName()
                clogp_sorted = sorted(clogPList)
                median_clogp = clogp_sorted[len(clogp_sorted)/2]
                min_clogp = clogp_sorted[0]
                max_clogp = clogp_sorted[len(clogp_sorted)-1]
                pki_sorted = sorted(pkiList)
                median_pki = pki_sorted[len(pki_sorted)/2]
                max_pki = pki_sorted[len(pki_sorted)-1]

                half = len(pki_sorted)/2
                adjustMax = (pki_sorted[len(pki_sorted)-1]-pki_sorted[(len(pki_sorted)+half)/2])/2
                if  adjustMax==0:
                    adjustMax = 0.1
                adjustmin = (pki_sorted[half/2]-pki_sorted[0])/2
                if adjustmin==0:
                    adjustmin = 0.1
                pki_axis = [pki_sorted[0]+adjustmin,pki_sorted[half/2],median_pki+0.1,pki_sorted[(len(pki_sorted)+half)/2]-0.1, pki_sorted[len(pki_sorted)-1]-adjustMax]
                collection = grob.VortexObjectCollection() #generate a graph collection
                x=None
                for k in range(vws.plotCount):
                    if isinstance(vws.getPlot(k),vp):
                        x = vws.getPlot(k)
                        break
                if x:
                    x.setXColumn(vtable.getColumnId(logp_column_name))
                    x.setYColumn(vtable.getColumnId("%s_%s(%s)"%(functionName,columnName,conc)))
                else:
                    vws.addChart(0,vtable.getColumnId(logp_column_name),vtable.getColumnId("%s_%s(%s)"%(functionName,columnName,conc)), -1)
                    x = vws.getPlot(0)
                max_x = max_clogp
                # min_x = x.getMinX()
                max_y = max_pki
                # min_y = x.getMinY()
                xLabel = 2
                i = 0
                for ic in range(2,7):
                    line = grob.VortexLine(-100,-100+ic,100,100+ic) #draw a line from x1,y1 to x2,y2
                    y_axis = max_y
                    x_axis = y_axis-ic
                    if x_axis > max_x:
                        x_axis = max_x
                        y_axis = x_axis+ic

                    # y_axis = pki_axis[i]
                    # x_axis = pki_axis[i]-ic
                    #
                    # if x_axis<min_clogp:
                    #     x_axis = min_clogp
                    #     y_axis = x_axis +ic
                    # if x_axis>max_clogp:
                    #     x_axis = max_clogp-2
                    #     y_axis=x_axis+ic
                    text = grob.VortexText('LiPE='+str(ic),x_axis-0.1,y_axis-0.1,0,0,0)
                    i+=1
                    text.setFillColor(Color.BLUE)
                    text.setFont(Font("Serif", Font.BOLD, 23))
                    collection.objects.add(text)
                    collection.objects.add(line) #add this line to the collection
               
                x.getBackgroundObjects().clear() #clear all background objs.
                x.getBackgroundObjects().add(collection) #add the collection to the plot
                x.dirtyRepaint() # repaint
                
            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setProgress(pair[0])
                self.progressMonitor.setNote(pair[1])


            def done(self):
                self.progressMonitor.close()
                try:
                    self.get()  #raise exception if abnormal completion
                    vtable.fireTableStructureChanged()
                    vortex.alert("Two new columns called %s, %s has been added to your Vortex table."%(self.newColumn.getName(),self.lipeCol.getName()))
                except ExecutionException, e:
                    vortex.alert(e.getMessage())

        task = ConvertTask()
        task.execute()
        
if __name__ == "__main__":
    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        lep = LEPlot()
        if lep.show_interface():
            lep.calculateAll()
