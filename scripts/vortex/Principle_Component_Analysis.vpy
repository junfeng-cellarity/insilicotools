#Imports of libraries
import os
import java.lang as lang
import java.lang.String as String
import java.net.URL as URL
import java.net.URLConnection as URLConnection
import java.io.InputStreamReader as InputStreamReader
import java.io.BufferedReader as BufferedReader
import java.io.DataOutputStream as DataOutputStream
#import javax.xml
import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.IntArray as IntArray
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.table.VortexJTable2 as VortexJTable2
import com.dotmatics.vortex.plots.VortexScatterPlot as VortexScatterPlot
import com.dotmatics.vortex.plots.VortexSymbol as VortexSymbol
import com.dotmatics.vortex.graphicalobjects as grob
import javax.swing.ProgressMonitor as ProgressMonitor
import javax.swing as swing
import java.awt.event.ActionListener as ActionListener
import javax.swing.event.DocumentListener as DocumentListener
import javax.swing.DefaultComboBoxModel as DefaultComboBoxModel
import java.lang.System as System
import java.util.Random as Random
import javax.swing.JScrollPane as JScrollPane
import java.awt.geom.GeneralPath as GeneralPath
import java.awt.geom.Line2D as Line2D
import java.awt.Color as Color
import java.awt.Font as Font
import java.lang.Math as Math
from javax.swing import BorderFactory, ProgressMonitor, SwingWorker, SwingUtilities
from java.lang import Runnable
from java.util.concurrent import ExecutionException

import sys

__doc__ = 'Principle Component Analysis using Existing Descriptors'
__author__ = 'Jun Feng'
__version__ = '1.0'

import com.dotmatics.vortex.util.Util as Util
import sys
sys.path.append(os.path.join(vortex.getVortexFolder(), "openeye"))
from initOEChem import VortexInit

vortexInit = VortexInit(vortex)
vortexInit.initOpeneye()
vortexInit.initJars()

from javastat.multivariate import PCA as PCA
import time
import java.net.URL
import math
import sys
import random


def isNanOrInf(value):
    if value == float('nan') or value == float('-inf') or value == float('inf'):
        return True
    return False

class PrincipleComponentAnalysis(lang.Runnable):
    def __init__(self):
        self.availableProperties = []
        self.selectedProperties = []
        self.selectedCompounds = []
        self.pcaLoadings = {}

        columns = vtable.getColumnsArray()
        for column in columns:
            self.availableProperties.append(column.getName())


    def show_interface(self):
        mypanel = swing.JPanel()

        mylabel = swing.JLabel("Select properties:")
        layout.fill(mypanel, mylabel, 0, 0)

        self.myduallist = components.VortexDualList()

        self.myduallist.setItems(self.availableProperties, self.selectedProperties)
        layout.fill(mypanel, self.myduallist, 0, 1)


        self.numPCACB = swing.JComboBox([2,3,4,5])
        tmpLbl2 = swing.JLabel("No. PCs:")
        subPanel = swing.JPanel(java.awt.FlowLayout(java.awt.FlowLayout.CENTER,20,10))
        subPanel.add(tmpLbl2)
        subPanel.add(self.numPCACB)

        layout.nofill(mypanel,subPanel,0,2)
        ret = vortex.showInDialog(mypanel, "Selection")

        if ret == vortex.OK:
            return True
        else:
            return False

    def calculateAll(self):
#        useSelected = self.mycheckbox.isSelected()
        self.useSelected = False
        self.nPCA = int(self.numPCACB.getSelectedItem())

        self.items = self.myduallist.getItems()
        if len(self.items) < 2:
            vortex.alert("At least 2 properties need to be selected.")
            return

        if len(self.items) < self.nPCA:
            vortex.alert("No more than %d components can be selected."%len(self.items))
            return

        for item in self.items:
            self.selectedProperties.append(self.availableProperties[item])

        class PrincipleComponentAnalysisTask(SwingWorker):
            def __init__(self,pca):
                self.progressMonitor = ProgressMonitor(vortex, "Progress","",0,100)
                self.progressMonitor.setMillisToDecideToPopup(0)
                self.progressMonitor.setMillisToPopup(0)
                self.pca = pca

            def process(self, pairs):
                pair = pairs[len(pairs)-1]
                self.progressMonitor.setProgress(pair[0])
                self.progressMonitor.setNote(pair[1])


            def done(self):
                self.progressMonitor.close()
                try:
                    percentExplained,newTable = self.get()  #raise exception if abnormal completion
                    vtable.fireTableStructureChanged()

                    propertyX = "PC_1"
                    propertyY = "PC_2"

                    column1 = newTable.findColumnWithName(propertyX,0)
                    column2 = newTable.findColumnWithName(propertyY,0)
                    column1id = newTable.getColumnId(column1)
                    column2id = newTable.getColumnId(column2)


                    vtable2 = VortexJTable2()
                    vtable2.setModel(newTable)
                    ws = vortex.createWorkspace(newTable,"%s_%s"%("PCA", vortex.getWorkspaceCount()))
                    ws.addChart(0,column1id,column2id)
                    ws.addComponent(JScrollPane(vtable2),"Table")
                    ws.setLayoutMode(VortexWorkspace.LAYOUT_MASTER)
                    plot = ws.getPlot(1)
                    plot.setDisplayHiddenPoints(True)

                    max_x = column1.getMax()+0.5
                    max_y = column2.getMax()+0.5

                    print >>sys.stderr,self.pca.pcaLoadings

                    xLoadings = self.pca.pcaLoadings[propertyX]
                    yLoadings = self.pca.pcaLoadings[propertyY]
                    graphCollections = grob.VortexObjectCollection()
                    for var in self.pca.selectedProperties:
                        print>> sys.stderr, var,xLoadings[var],yLoadings[var],max_x,max_y
                        self.drawLoading(graphCollections,var,xLoadings[var],yLoadings[var],max_x,max_y)
                    plot.getBackgroundObjects().clear()
                    plot.getBackgroundObjects().add(graphCollections)
                    plot.dirtyRepaint()
                    ws.updateLayout(VortexWorkspace.LAYOUT_MASTER)

                    class myRunnalbe(Runnable):
                        def run(self):
                            plot.setMaxX(column1.getMax()+2.0)
                            plot.setMaxY(column2.getMax()+2.0)

                    SwingUtilities.invokeLater(myRunnalbe())
                    vortex.alert("Percent Explained:%5.1f"%percentExplained)

                except ExecutionException, e:
                    vortex.alert(e.getMessage())

            def drawLoading(self,collection,varName, xloadings, yloadings,max_x, max_y):
                # add a line x1,y1, x2,y2
                offset = (max_x,max_y)
                line = grob.VortexLine(offset[0],offset[1], xloadings+offset[0],yloadings+offset[1])
                collection.objects.add(line)

                startx = offset[0]
                starty = offset[1]
                endx = xloadings+offset[0]
                endy = yloadings+offset[1]

                arrowLength = Math.sqrt((startx-endx)*(startx-endx)+(starty-endy)*(starty-endy))/8.0
                arrowWidth = arrowLength/3.0
                angle = Math.atan((endy-starty)/(endx-startx))

                sigma = -1 if endx>startx else 1

                arrowBaseX = endx + sigma*Math.cos(angle) * arrowLength
                arrowBaseY = endy + sigma*Math.sin(angle) * arrowLength

                labelX = endx - sigma*Math.cos(angle) * arrowLength
                labelY = endy - sigma*Math.sin(angle) * arrowLength

                arrowLeftX = arrowBaseX+ Math.cos(angle + Math.PI / 2.0) * arrowWidth
                arrowLeftY = arrowBaseY+ Math.sin(angle + Math.PI / 2.0) * arrowWidth

                arrowRightX = arrowBaseX- Math.cos(angle + Math.PI / 2.0) * arrowWidth
                arrowRightY = arrowBaseY- Math.sin(angle + Math.PI / 2.0) * arrowWidth

                line = grob.VortexLine(endx,endy,arrowLeftX,arrowLeftY)
                collection.objects.add(line)

                line = grob.VortexLine(endx,endy,arrowRightX,arrowRightY)
                collection.objects.add(line)



                #add it to the collection

                # create a text object
                # text value, x, y, rotatation, xfrac offset, yfrac offset
                text = grob.VortexText(varName, labelX,labelY, 0, 0, 0)

                # set the text fill color
                text.setFillColor(Color.RED)

                # set the text font
                text.setFont(Font("Serif", Font.PLAIN, 8))

                # add it to the collection
                collection.objects.add(text)


            def doInBackground(self):

                nrow = vtable.getRealRowCount()
                rows = range(0,nrow)
                if self.pca.useSelected:
                    nrow = vtable.getSelectedRowCount()
                    rows = []
                    tmp = vtable.getSelectedRows()
                    for i in range(0,tmp.size()):
                        rows.append(tmp.get(i))

                data = []
                columns = []

                self.super__publish([(1,"Normalizing input data ...")])

                for colName in self.pca.selectedProperties:
                    column = vtable.findColumnWithName(colName,0)
                    if 0 != column.getStdDev():
                        columns.append(column)
                        rowData = []
                        for row in range(0,nrow):
                            rawValue = float(column.getValue(rows[row]))
                            if isNanOrInf(rawValue):
                                rawValue = 0
                            rowData.append(math.fabs((rawValue-column.getMean())/column.getStdDev()))
                        data.append(rowData)

                self.super__publish([(5,"Calculating Principle Components ...")])
                #Doing Principle component analysis with 100% (1.0) explained
                pca = PCA(1.0)
                principle_components = pca.principalComponents(data)
                eigenvalues = pca.getEigenvalues()
                eigenvalues.reverse()

                self.super__publish([(25,"Building new columns ...")])
#                for columnId in vtable.getAllColumnIds():
#                    column = vtable.getColumn(columnId)
#                    if(column.getName().startswith("PC")):
#                        for i in range(0,vtable.getRealRowCount()):
#                            column.setType(Column.COLUMN_DOUBLE)
#                            column.setDouble(i,0.0)

                if self.pca.nPCA > len(principle_components):
                    self.pca.nPCA = len(principle_components)

                sum = 0.0
                sum1 = 0.0
                for idx,value in enumerate(eigenvalues):
                    sum += value
                    if idx < self.pca.nPCA:
                        sum1 += value
                percentExplained = 100.0*sum1/sum
                # print eigenvalues
                # print "Percent explained:%5.1f"%percentExplained

                for k in range(0,self.pca.nPCA):
                    pcaColumnName = "PC_%d"%(k+1)
                    self.pca.pcaLoadings[pcaColumnName] = {}
                    for id,colName in enumerate(self.pca.selectedProperties):
                        self.pca.pcaLoadings[pcaColumnName][colName] = principle_components[k][id]

                new_columns = []
                newTable = VortexTableModel()
                for column in vtable.getColumnsArray():
                    if column.getName()=="Structure":
                        new_column = newTable.findColumnWithName("MOLFILE",1,Column.COLUMN_STRING)
                    else:
                        new_column = newTable.findColumnWithName(column.getName(),1, column.getType())
                    for i in range(0,nrow):
                        if column.getName()=="Structure":
                            molFile = vtable.getMolFileManager().getMolFileAtRow(rows[i])
                            new_column.setString(rows[i],molFile)
                        else:
                            value = column.getValueAsString(rows[i])
                            if value!="N/A":
                                new_column.setValueFromString(i,value)
                            else:
                                new_column.setValueFromString(i,None)

                    new_columns.append(new_column)

                for k in range(0,self.pca.nPCA):
                    pcaColumnName = "PC_%d"%(k+1)
                    new_column = newTable.findColumnWithName(pcaColumnName,1,Column.COLUMN_DOUBLE)
                    for i in range(0,nrow):
                        pcaValue = 0.0
                        for j in range(0,len(principle_components)):
                            pcaValue += principle_components[k][j]*data[j][i]
                        new_column.setDouble(i,pcaValue)
                return percentExplained,newTable

        task = PrincipleComponentAnalysisTask(self)
        task.execute()

if __name__ == "__main__":

    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        pca = PrincipleComponentAnalysis()
        if pca.show_interface():
            pca.calculateAll()