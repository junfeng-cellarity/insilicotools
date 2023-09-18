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
import javax.swing.JComboBox as JComboBox
import javax.swing.JTextField as JTextField
import javax.swing.JButton as JButton
import javax.swing.JDialog as JDialog
import java.awt.event.ActionListener as ActionListener
import javax.swing.event.DocumentListener as DocumentListener
import javax.swing.JLabel as JLabel
import javax.swing.JPanel as JPanel
import javax.swing.SwingUtilities as SwingUtilities
import javax.swing.BorderFactory as BorderFactory
import java.awt.BorderLayout as BorderLayout
import java.awt.GridLayout as GridLayout
import java.awt.Dimension as Dimension
import javax.swing.JScrollPane as JScrollPane
import javax.swing.JSplitPane as JSplitPane
import javax.swing.JPanel as JPanel

import com.dotmatics.vortex.components as components
import com.dotmatics.vortex.table.Column as Column
import com.dotmatics.vortex.util.Util as Util
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.IntArray as IntArray
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.plots.VortexScatterPlot as VortexScatterPlot

__doc__ = 'Add Compound Tracking Columns'
__author__ = 'Jun Feng'
__version__ = '1.0'

import datetime

columnPickListDict = {
    "Vendor":"Chempartner,Pharmaron,Aragen",
    "Status":"In Progress,On Hold,Redesign,Selected,Complete",
    "Priority":"N/A,High,Low,Medium,Top",
    "Contact":"Brian,Mike,Tiansheng,Jotham,O'Neill,Siva,Dan,Ingo",
    "Assignment":"N/A,Team 1,Team 2,Team 3,Team 4,Team 5,Team 6,Team 7,Team 8,Team A, Team C, Team S, Team FP, Unassigned,Flex,3371,4269,0978,7349,2825"
}


if __name__ == "__main__":
    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        mypanel = JPanel()
        mypanel.setLayout(GridLayout(len(columnPickListDict.keys()),1))
        comboList = {}
        for key in columnPickListDict.keys():
            cb = JComboBox(columnPickListDict[key].split(","))
            label = JLabel("%s:"%key)
            p = JPanel()
            p.add(label)
            p.add(cb)
            mypanel.add(p)
            comboList[key] = cb
        ret = vortex.showInDialog(JScrollPane(mypanel), "Pick default value for each column:")
        if ret == vortex.OK:
            date_col = vtable.findColumnWithName("Date",True,Column.COLUMN_STRING)
            for colname in columnPickListDict.keys():
                default_value = comboList[colname].getSelectedItem()
                vtable.addPickListColumn(default_value,columnPickListDict[colname],True,colname)
                col = vtable.findColumnWithName(colname,False)
                for row in range(0,vtable.getRealRowCount()):
                    col.setValueFromString(row, default_value)
            date = datetime.date.today().isoformat()
            for row in range(0,vtable.getRealRowCount()):
                date_col.setValueFromString(row,date)
            vtable.fireTableStructureChanged()



