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
import com.dotmatics.vortex.Vortex as Vortex
import com.dotmatics.vortex.table.IntArray as IntArray
import com.dotmatics.vortex.VortexWorkspace as VortexWorkspace
import com.dotmatics.vortex.plots.VortexScatterPlot as VortexScatterPlot

__doc__ = 'Add SMARTS Pattern match count as a calculated column'
__author__ = 'Jun Feng'
__version__ = '1.0'

import sys
import time
import java.net.URL
import math
import jarray
import os

from java.lang import Runnable, Double
from java.util import Random, BitSet
from java.util.concurrent import ExecutionException

from javax.swing.border import TitledBorder
from javax.swing import BorderFactory, ProgressMonitor, JTextArea, JPanel, JLabel, JCheckBox
from javax.swing import JButton, JFrame, JTextField,JComboBox
from java.awt import GridBagLayout as awtGridBagLayout, Dimension, BorderLayout
from java.awt import GridBagConstraints as awtGridBagConstraints
from java.awt import Insets as awtInsets
import javax.swing.SwingWorker as SwingWorker



if __name__ == "__main__":
    if vws is None:
        vortex.alert("You must have a workspace loaded...")
    else:
        dict = {}
        col1 = vtable.findColumnWithName("CLUSTER",0,Column.COLUMN_INT)
        col2 = vtable.findColumnWithName("DISTANCE",0,Column.COLUMN_DOUBLE)
        nrow = vtable.getRealRowCount()
        rows = range(0,nrow)
        for row in rows:
            vtable.setRowSelected(row,False)
            cluster = col1.getInt(row)
            distance = col2.getDouble(row)
            if cluster not in dict and distance == 0:
                vtable.setRowSelected(row,True)
                dict[cluster]=1
        vtable.fireTableStructureChanged()




