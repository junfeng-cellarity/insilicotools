# import openpyxl
# wb = openpyxl.load_workbook('node_interactions.xlsx')
# sheet = wb.get_sheet_by_name('node_interactions')

from Tkinter import *
import Tkinter
canvas = Canvas(width=800, height= 800, bg='white')  
canvas.pack(expand=YES, fill=BOTH)

protiens = []
protiens2 = []
num = 0
num1 = 0
start_x = 500 #initial position
start_y = 200
#num2 = 0

# for cell in sheet.columns[0]:
#     protiens.append(cell.value)
#
# for cell in sheet.columns[1]:
#     protiens2.append(cell.value)

# protiens.remove('node1')
# protiens2.remove('node2')

class Node():
    def __init__(self, target):
        self.target = target
        self.parent = None
        self.children = []
        self.x = start_x
        self.y = start_y

    def numOfChildren(self):
        print("Total number of Nodes are: ", len(self.children))

    def calculateX(self,idx,total,x1,x2):
        div = (x2-x1)/total
        print div
        return x1+idx*div+div/2

    def addChild(self, childNode): # recalculate child x coordinates, when new child added
        self.children.append(childNode)
        childNode.parent = self
        childNode.y = 100+self.y
        total_span = 100*len(self.children)
        x1 = self.x - total_span/2
        x2 = self.x + total_span/2
        for id,child in enumerate(self.children):
            child.x = self.calculateX(id,len(self.children),x1,x2)


    def displayNode(self):
        x1 = self.x-30
        y1 = self.y-20
        x2 = self.x+30
        y2 = self.y+20
        canvas.create_oval(x1,y1,x2,y2,fill="blue")
        for child in self.children:
            canvas.create_line(self.x,self.y,child.x,child.y,fill="black",arrow=LAST)
            child.displayNode()


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
