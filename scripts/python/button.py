import Tkinter as tk
root = tk.Tk()

def myfunction(event):
    print buttons[event.widget]

buttons = {}
for i in range(10):
    b = tk.Button(root, text='button' + str(i))
    buttons[b] = i # save button, index as key-value pair
    b.bind("<Button-1>", myfunction)
    b.place(x=10,y=(10+(25*i)))
root.mainloop()

