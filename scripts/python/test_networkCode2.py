"""import the package to manipulate the cells directly"""

import openpyxl
wb = openpyxl.load_workbook('node_interactions3.xlsx')
sheet = wb.get_sheet_by_name('node_interactions')

"""import the graphics package to print the circles"""
from Tkinter import *
canvas = Canvas(width=800, height= 800, bg='white')
canvas.pack(expand=YES, fill=BOTH)

"""importing sorting library"""
from collections import OrderedDict

"""iterator and temporary variables"""
i = 0
child = ""
existingNodes = {}


"""initialize the lists that contain column 1 and column2's contents"""
targets1 = []
targets2 = []

start_x = 400
#start_y = 200
start_y = 20

"""put all the values from the first column into targets1"""
for cell in sheet.columns[0]:
    targets1.append(str(cell.value))

"""put all the values from the second column into targets2"""
for cell in sheet.columns[1]:
    targets2.append(str(cell.value))

"""remove node1 and node2"""
targets1.remove('node1')
targets2.remove('node2')

print(targets1)
print(targets2)

maxSpan = 100

class Node():
    def __init__(self, target):
        self.target = target
        self.parent = None
        self.children = []
        self.x = start_x
        self.y = start_y
        self.layer = 0
        self.visible = False

    def numOfChildren(self):
        print("Total number of Nodes are: ", len(self.children))

    def calculateX(self,idx,total,x1,x2):
        div = (x2-x1)/total
        print (div)
        return x1+idx*div+div/2

    def addChild(self, childNode): # recalculate child x coordinates, when new child added
        self.children.append(childNode)
        childNode.parent = self
        if childNode.layer == 0:
            childNode.layer = self.layer +1
        if childNode.y==start_y:
            childNode.y = 100+self.y
        total_span = 100*len(self.children)
        x1 = self.x - total_span/2
        x2 = self.x + total_span/2
        for id,child in enumerate(self.children):
            if child.x==start_x:
                child.x = self.calculateX(id,len(self.children),x1,x2)


    def displayNode(self):
        if self.visible:
            return
        self.visible = True
        #total_span = 100*len(self.children)
        #x_1 = self.x - total_span/2
        #x_2 = self.x + total_span/2
        #for id,child in enumerate(self.children):
#            child.x = self.calculateX(id,len(self.children),x_1,x_2)
        x1 = self.x-30
        y1 = self.y-20
        x2 = self.x+30
        y2 = self.y+20
        canvas.create_oval(x1,y1,x2,y2,fill="light blue")
        canvas.create_text(self.x,self.y, text = self.target)
        for child in self.children:
            canvas.create_line(self.x,self.y,child.x,child.y,fill="black",arrow=LAST)
            child.displayNode()
            #print(len(self.children))

target_list = []
for target in targets1:
    if target not in target_list:
        target_list.append(target)
        existingNodes[target] = Node(target)

for target in targets2:
    if target not in target_list:
        target_list.append(target)
        existingNodes[target] = Node(target)

print target_list
print existingNodes

for id,target in enumerate(targets1):
    existingNodes[target].addChild(existingNodes[targets2[id]])



existingNodes[targets1[0]].displayNode()
top = Tk()
top.mainloop()


"""        
root = Node("root")

child1= Node("child1")
root.addChild(child1)

child2 = Node("grandchild1")
child1.addChild(child2)

child3 = Node("grandgrandchild1")
child2.addChild(child3)
child4 = Node("grandgrandchild2")
child2.addChild(child4)

child5 = Node("gggchild1")
child6 = Node("gggchild1")
child7 = Node("gggchild1")
child8 = Node("gggchild1")
child3.addChild(child5)
child3.addChild(child6)
child3.addChild(child7)
child3.addChild(child8)

root.displayNode()
top = Tkinter.Tk()
top.mainloop()
"""
