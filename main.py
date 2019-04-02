from imports import *
from utils import *
from BackEndGame import *
from FrontEndGame import *

pa = pickle.load(open('Resources/Network_pypath_burn.by', 'br'))
root=tk.Tk()
gui=Gui(root, pa)
root.mainloop()