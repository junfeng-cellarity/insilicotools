import openpyxl
wb = openpyxl.load_workbook('node_interactions.xlsx')
sheet = wb.get_sheet_by_name('node_interactions')

from Tkinter import *

canvas = Canvas(width=800, height= 600, bg='white')  
canvas.pack(expand=YES, fill=BOTH)

num = 0
num1 = 0
num2 =50
i = 0
celVal = ''
celVal2 = ''
exists = ''
protiens = []
newProtiens = []
protiens2 = []

for cell in sheet.columns[0]:
    canvas.create_oval(250, (5+num), 300, (55+num), fill="light blue")
    #canvas.create_oval(5, (5+num), 50, (50+num), fill="blue")
    canvas.create_text((275),(5+22+num),text = cell.value)
    num += 55
    i += 1

#put all the values in the first column into protiens
for cell in sheet.columns[0]:
    protiens.append(cell.value)

#put all the values in the second column into protiens2
for cell in sheet.columns[1]:
    protiens2.append(cell.value)

#delete all the values in protiens 2 that are in protiens
"""for n in protiens:
    if n in protiens2:
        protiens2.remove(n)
"""
#remove node1 and node2
protiens.remove('node1')
protiens2.remove('node2')

#delete repeated values in protiens
for n in protiens:
    if n not in newProtiens:
        newProtiens.append(n)

    

#print(newProtiens)
#print(protiens2)

#turn protiens into a series of bubbles with the protien names in the bubbles
for n in newProtiens:
    canvas.create_oval(250, (5+num1), 300, (55+num1), fill="light blue")
    canvas.create_line(275, 55+num1, 275, 67+num1, fill = "black", arrow = LAST) #creates a line between the bubbles
    canvas.create_text((275),(5+25+num1), text = n)

    if (protiens2[i] != protiens[i+1]):
        canvas.create_oval((250 + num2), (5+num1), (300+num2), 55+num1, fill="light blue")
        #canvas.create_text((275),(5+25+num), text = protiens2[i])
 
    num1 += 60
    i+=1
    

canvas.create_oval(250, (5+num1), 300, (55+num1), fill="light blue")
canvas.create_text((275),(5+22+num1), text = sheet['B13'].value)


#manipulating cells
for cell in sheet.columns[0]:
    if (cell.value == 'node1'):
        num1 -= 52.3
    if (cell.value != celVal):
        canvas.create_oval(250, (5+num1), 300, (55+num1), fill="light blue")
        canvas.create_text((275),(5+22+num1),text = cell.value)
    else:
        num1 -= 50
        
    celVal = cell.value
    num1 += 50
    i += 1
    
canvas.create_oval(250, (5+num1), 300, (55+num1), fill="light blue")
canvas.create_text((275),(5+22+num1),text = sheet['B' + str(i)].value)



for cell in sheet.columns[1]:
    if (cell.value == 'node2'):
        num2 -= 102
    if (cell.value != celVal2):
        canvas.create_oval(400, (5+num2), 350, (55+num2), fill="light blue")
        canvas.create_text((375),(5+22+num2),text = cell.value)
    else:
        num2 -= 50
    celVal2 = cell.value
    num2 += 50
    i += 1


        
    
        
