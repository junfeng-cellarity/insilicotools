import sys
import com.dotmatics.vortex.util.Util as Util


import collections
import org.apache.commons.math.stat.correlation.PearsonsCorrelation as pc
import com.dotmatics.vortex.plots.VortexScatterPlot as vp
import com.dotmatics.vortex.components as comp
from javax.swing import JPanel,BoxLayout,JScrollPane,JFrame,SwingUtilities, JTable, UIManager,SwingWorker, ProgressMonitor,JTextField,RowFilter,SortOrder
from javax.swing.table import AbstractTableModel, TableRowSorter
from javax.swing.event import ListSelectionEvent,ListSelectionListener,DocumentListener
from java.awt.event import ActionEvent,ActionListener,MouseAdapter,MouseListener
from java.awt import GridLayout,Dimension
from java.lang import Runnable
from java.util import Comparator

class DualListDialog:

    def __init__(self):
        self.success = False
        self.dualList = None
        self._buildDualList()

    def _buildDualList(self):
        leftList = []
        rightList = []
        self.dualList = comp.VortexDualList()
        self.dualList.setSourceString('Columns')
        self.dualList.setDestString('Selected Columns')
        forChoose = []
        self.colDict = {}
        ct = 0
        for idx in xrange(vtable.getColumnCount()):
            col = vtable.getColumn(idx)
            if col.getForceNumeric() == False and col.isNumeric() != True:
                continue
            self.colDict[ct] = idx
            forChoose.append(col.getNiceName())
            ct += 1
        for colName in forChoose:
            leftList.append(colName)
        self.dualList.setItems(leftList, rightList)

    def show(self):
        resp = vortex.showInDialog(self.dualList, "Correlation Column Selection")
        if resp == vortex.OK:
            if len(self.dualList.getItems())<2:
                vortex.alert('Please choose at least two columns')
                self.show()
            else:
                self.success = True


class MyTableModel (AbstractTableModel):
    def __init__(self, data=None):
        self.columnNames = ("Col1",  "Col2", "R","R2","Size")
        self.data = data
        if not data:
            self.data = []

    def setData(self, data):
        self.data = data

    def getColumnCount(self):
        return len(self.columnNames)

    def getRowCount(self):
        return len(self.data)

    def getColumnName(self, col):
        return self.columnNames[col]

    def getValueAt(self, row, col):
        return self.data[row][col]


class myComparator(Comparator):

    def compare(self, o1, o2):
        so1 = str(float(o1))
        so2 = str(float(o2))
        if so1 == 'nan' or so2 == 'nan':
            if o1 < o2:
                return -1
            elif o1 == o2:
                return 0
            else:
                return 1
        else:
            if float(o1) < float(o2):
                return -1
            elif float(o1) == float(o2):
                return 0
            else:
                return 1


class CorrelationTable(JPanel):
    def __init__(self):
        global data
        try:
            self.setLayout(BoxLayout(self,BoxLayout.Y_AXIS))
            self.table = JTable(MyTableModel(data))
            self.table.setPreferredScrollableViewportSize(Dimension(500, 500));
            self.table.setFillsViewportHeight(True);
            self.table.getSelectionModel().addListSelectionListener(self.RowListener(self))
            self.rowSorter = TableRowSorter(self.table.getModel())
            for i in range(2,5):
                self.rowSorter.setComparator(i,myComparator())
            self.table.setRowSorter(self.rowSorter)
            self.table.getTableHeader().addMouseListener(self.ColumnHeaderListener(self))
            self.filterTextBox = JTextField('',30)
            self.filterTextBox.getDocument().addDocumentListener(self.myDocumentListener(self))
            self.add(self.filterTextBox)
            self.add(JScrollPane(self.table))
            self.prevPlot=[]
        except Exception,e:
            vortex.alert(str(e))

    def textFilter(self):
        rf = None
        try:
            rf = RowFilter.regexFilter('(?i)'+self.filterTextBox.getText(),[0,1]) #(?i) to make it case insensitive
            self.rowSorter.setRowFilter(rf)
        except Exception,e:
            vortex.alert(str(e))

    class ColumnHeaderListener(MouseAdapter):
        def __init__(self,parent):
            self.parent=parent
        def mouseClicked(self,evt):
            try:
                colModel = self.parent.table.getColumnModel()
                viewColIndex = colModel.getColumnIndexAtX(evt.getX())
                modelColIndex = self.parent.table.convertColumnIndexToModel(viewColIndex)
                if (viewColIndex==-1):
                    return
                if modelColIndex==0:
                    if len(self.parent.table.getRowSorter().getSortKeys())==2 and self.parent.table.getRowSorter().getSortKeys()[0].getColumn()==0:
                        return
                    self.parent.table.getRowSorter().setSortKeys([TableRowSorter.SortKey(0,SortOrder.ASCENDING),TableRowSorter.SortKey(3,SortOrder.DESCENDING)])

                elif modelColIndex==1:
                    if len(self.parent.table.getRowSorter().getSortKeys())==2 and self.parent.table.getRowSorter().getSortKeys()[0].getColumn()==1:
                        return
                    self.parent.table.getRowSorter().setSortKeys([TableRowSorter.SortKey(1,SortOrder.ASCENDING),TableRowSorter.SortKey(3,SortOrder.DESCENDING)])
                elif modelColIndex==2:
                    if len(self.parent.table.getRowSorter().getSortKeys())==3 and self.parent.table.getRowSorter().getSortKeys()[0].getColumn()==2:
                        return
                    self.parent.table.getRowSorter().setSortKeys([TableRowSorter.SortKey(2,SortOrder.DESCENDING),TableRowSorter.SortKey(0,SortOrder.ASCENDING),TableRowSorter.SortKey(1,SortOrder.ASCENDING)])
                elif modelColIndex==3:
                    if len(self.parent.table.getRowSorter().getSortKeys())==3 and self.parent.table.getRowSorter().getSortKeys()[0].getColumn()==3:
                        return
                    self.parent.table.getRowSorter().setSortKeys([TableRowSorter.SortKey(3,SortOrder.DESCENDING),TableRowSorter.SortKey(0,SortOrder.ASCENDING),TableRowSorter.SortKey(1,SortOrder.ASCENDING)])
                elif modelColIndex==4:
                    if len(self.parent.table.getRowSorter().getSortKeys())==3 and self.parent.table.getRowSorter().getSortKeys()[0].getColumn()==4:
                        return
                    self.parent.table.getRowSorter().setSortKeys([TableRowSorter.SortKey(4,SortOrder.DESCENDING),TableRowSorter.SortKey(3,SortOrder.DESCENDING),TableRowSorter.SortKey(0,SortOrder.ASCENDING),TableRowSorter.SortKey(1,SortOrder.ASCENDING)])
            except Exception, e:
                vortex.alert(str(e))

    class myDocumentListener(DocumentListener):
        def __init__(self,parent):
            self.parent = parent
        def changedUpdate(self,event):
            try:
                self.parent.textFilter()
            except Exception, e:
                vortex.alert(str(e))
        def insertUpdate(self,event):
            try:
                self.parent.textFilter()
            except Exception, e:
                vortex.alert(str(e))
        def removeUpdate(self,event):
            try:
                self.parent.textFilter()
            except Exception, e:
                vortex.alert(str(e))

    class RowListener (ListSelectionListener):
        def __init__(self,parent):
            self.parent =parent
        def valueChanged(self,event):
            if event.getValueIsAdjusting():
                return
            try:
                if len(self.parent.prevPlot):
                    for plot in self.parent.prevPlot:
                        plot.dispose()
                for row in self.parent.table.getSelectedRows():
                    convertedRow = self.parent.table.convertRowIndexToModel(row)
                    colI = self.parent.table.getModel().getValueAt(convertedRow,0)
                    colJ = self.parent.table.getModel().getValueAt(convertedRow,1)

                    idxI = self.parent.table.convertColumnIndexToView(0)
                    idxJ = self.parent.table.convertColumnIndexToView(1)
                    print str(idxI>idxJ)
                    x=None
                    for k in range(vws.plotCount):
                        if isinstance(vws.getPlot(k),vp):
                            x = vws.getPlot(k)
                            break
                    if x:
                        if idxI>idxJ:
                            x.setXColumn(vtable.getColumnId(colI))
                            x.setYColumn(vtable.getColumnId(colJ))
                        else:
                            x.setXColumn(vtable.getColumnId(colJ))
                            x.setYColumn(vtable.getColumnId(colI))
                    else:
                        if idxI>idxJ:
                            chart = vws.addChart(0,vtable.getColumnId(colI),vtable.getColumnId(colJ), -1)
                        else:
                            chart = vws.addChart(0,vtable.getColumnId(colJ),vtable.getColumnId(colI), -1)
                        self.parent.prevPlot.append(chart)
            except Exception,e:
                print e


class CorrelationTask(SwingWorker):
    def __init__(self, data, dualList, colDict):
        self.colDict = colDict
        self.dualList = dualList
        self.progressMonitor = ProgressMonitor(vortex, "Calculating Correlation ...", "", 0, 100)
        self.data = data
        self.succeeded = True

    def doInBackground(self):
        p = pc()
        rowCount = vtable.getRealRowCount()
        data = []
        dataDict = collections.defaultdict(lambda: collections.defaultdict(lambda: None))
        self.super__publish([(0/vtable.getColumnCount(), "process column %d" % (0,))])
        items = self.dualList.getItems()
        for idx in xrange(len(items)):
            col = vtable.getColumn(self.colDict[items[idx]])
            m = 0
            while m < rowCount:
                vI = col.getValueAsString(m)
                if vI is None or vI == 'N/A':
                    m += 1
                    continue
                dataDict[items[idx]][m] = col.getValue(m)
                m += 1
        print 'here'
        for i in xrange(len(items)):
            if self.progressMonitor.isCanceled():
                break
            if items[i] not in dataDict:
                print items[i],"not found"
                continue
            self.super__publish([(100*i/len(items),"process column %d"%(i))])
            colI = vtable.getColumn(self.colDict[items[i]])
            colNameI = colI.getNiceName()
            for j in xrange(i+1,len(items)):
                if items[j] not in dataDict:
                    continue
                colJ = vtable.getColumn(self.colDict[items[j]])
                colNameJ = colJ.getNiceName()
                n=0
                listINew =[]
                listJNew =[]
                for idx in range(rowCount):
                    if dataDict[items[i]][idx] is None or dataDict[items[j]][idx] is None:
                        continue
                    listINew.append(dataDict[items[i]][idx])
                    listJNew.append(dataDict[items[j]][idx])

                if len(listINew)<2 or len(listJNew)<2:
                    continue
                R  = p.correlation(listINew,listJNew)
                if str(R)=='nan':
                    continue
                R2 = R**2
                data.append([colNameI,colNameJ,round(R,2),round(R2,2),min(len(listINew),len(listJNew))])


        self.data.extend(sorted(data,key = lambda row:row[3],reverse=True))

    def process(self, pairs):
        pair = pairs[len(pairs)-1]
        self.progressMonitor.setProgress(pair[0])
        self.progressMonitor.setNote(pair[1])

    def done(self):
        self.progressMonitor.close()
        #after done, show the Correlation Table
        if self.succeeded:
            newContentPane = CorrelationTable()
            vortex.showOKDialog(newContentPane, "Correlation Table", None)
            try:
                self.get()  #raise exception if abnormal completion
            except Exception, e:
                vortex.alert(e.getMessage())

if __name__ == "__main__":
    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        data = []
        dualListDialog = DualListDialog()
        dualListDialog.show()
        if dualListDialog.success:
            task = CorrelationTask(data, dualListDialog.dualList, dualListDialog.colDict)
            task.execute()

