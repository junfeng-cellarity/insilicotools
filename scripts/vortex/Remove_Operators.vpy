__author__ = 'jfeng1'

import sys
import com.dotmatics.vortex.table.VortexTableModel as VortexTableModel
import com.dotmatics.vortex.util.Layout as layout
import com.dotmatics.vortex.components as components

columnNames = vtable.getColumnNames(0)
cols = []
includeCols = []

for name in columnNames:
    cols.append(name)

# Prompt user for list of inhibition values to de-pivot
content = javax.swing.JPanel()
duallist = components.VortexDualList()
duallist.setItems(cols, None)
layout.fill(content, duallist, 0, 0)
ret = vortex.showInDialog(content, "Choose the columns to fix:")

if ret == vortex.OK:
    for i in duallist.getItems():
        includeCols.append(cols[i])

for col in includeCols:
    cols.remove(col)

if not len(includeCols):
    vortex.alert("No columns selected")
    sys.exit(0)

#vtable = VortexTableModel()
for col in includeCols:
    column = vtable.findColumnWithName(col,0)
    for i in range(0,vtable.getRowCount()):
        v = column.getValueAsString(i)
        if v is not None:
            v = v.replace('GT','').replace('>','').replace('<','')
            column.setValueFromString(i,v)
    column.setForceNumeric(True)

vortex.alert("finished.")