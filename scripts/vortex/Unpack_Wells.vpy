# Unpack well values from a plate dataset
# Well values are assumed to be of the form letters followed by digits
# A 'Row' column is created for the letter part
# and a 'Column' column is created for the digit part
#
# Author: Mike Hartshorn
# Copyright (C) dotmatics Limited, 2008

#import re
import java.util.regex.Matcher as Matcher
import java.util.regex.Pattern as Pattern

rowlabel 	= swing.JLabel("Row name")
collabel 	= swing.JLabel("Column name")
rowtf 		= swing.JTextField("Row")
coltf 		= swing.JTextField("Column")
welllabel 	= swing.JLabel("Well column")
wellcb 		= workspace.getColumnComboBox()

panel 		= swing.JPanel()

layout.fill(panel, welllabel, 	0, 0)
layout.fill(panel, wellcb, 	1, 0)
layout.fill(panel, rowlabel, 	0, 1)
layout.fill(panel, rowtf, 	1, 1)
layout.fill(panel, collabel, 	0, 2)
layout.fill(panel, coltf, 	1, 2)

ret = vortex.showInDialog(panel, "Unpack Well Values")

if ret == vortex.OK:
	wellidx = wellcb.getSelectedIndex()
	
	if wellidx == 0:
		vortex.alert("you must choose a column")
	else:
		wellcol = vtable.getColumn(wellidx - 1)
	
		rowname = rowtf.getText()
		colname = coltf.getText()
		
		rowcol = vtable.findColumnWithName(rowname, 1)
		rowcol.setType(vortex.STRING)
		colcol = vtable.findColumnWithName(colname, 1)
		colcol.setType(vortex.INT)
		
		rows = vtable.getRealRowCount()

		#p = re.compile("([A-Za-z]+)([0-9]+)")
		p = Pattern.compile("([A-Za-z]+)([0-9]+)")
		for r in range(0,rows):
			wellval = wellcol.getValueAsString(r)
			
			#m = p.match(wellval)
			m=p.matcher(wellval)
			m.matches()

			#print wellval, m.group(1), m.group(2)
			
			rowcol.setValueFromString(r, m.group(1))
			colcol.setValueFromString(r, m.group(2))
	
	vtable.fireTableStructureChanged()