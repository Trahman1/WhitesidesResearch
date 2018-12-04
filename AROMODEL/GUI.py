from Tkinter import *
import Configure
import os
import ovitoHandler
import Logger
import Controller
from PIL import ImageTk
from PIL import Image
import GUI.GUIHelper

window = Tk()
window.title("Aromodel")
RESOLUTION = '1000x600'

class listBoxFrame(Frame):
	def __init__(self, master, title, **keywords):
		Frame.__init__(self, master, relief=SUNKEN, width=400, height=200, padx=25)
		# self.grid(column = 0, row = 0)
		Label(self, text=title, font = ("Arial Bold", 20)).grid(column=0, row=0)
		self.listbox = Listbox(self, height=8)#, yscrollcommand = statusscroll.set)#
		self.listbox.grid(column=0, row=2, sticky = N+S+E+W)
		# dict takes input to output value
		self.dict = {}
		self.dict[""]=None

	def insert(self, i,text,selection_output):
		self.listbox.insert(i,text)
		self.dict[text]=selection_output

	def setOnselect(self, OnSelect):
		__onselect = lambda self=self: OnSelect()
		self.listbox.bind('<<ListboxSelect>>', __onselect)

	def selection(self):
		return self.dict[self.listbox.get(ANCHOR)]

	def clear(self):
		self.listbox.delete(0, END)

	def refreshList(self, lst):
		self.clear()
		for text, output in lst:
			self.insert(END,text,output)

class plotFrame(Frame):
	def __init__(self, master, **keywords):
		Frame.__init__(self, master, relief=SUNKEN, bg = "gray", width=500, height=600, padx=25)
		self.grid_propagate(False)
		Label(self, text="Plots", bg="gray", font = ("Arial Bold", 20), padx=25).grid(column=0, row=0, sticky=N)
		Button(self, text="Tangent Correlation",
			command = lambda frame=self: frame.addTangentCorrelation()).grid(column=0, row=1)
		Button(self, text="Carbon Density",
			command = lambda frame=self: frame.addCarbonDensity()).grid(column=0, row=2)
		Button(self, text="Pair Correlation",
			command = lambda frame=self: frame.addPairCorrelation()).grid(column=0, row=3)

		self.plots = [None]
		self.plotresx, self.plotresy = 400., 400.
		self.plotlabels = []
		self.plotimages = []

	def setOnselect(self, OnSelect):
		return lambda self=self: OnSelect()

	def refresh(self):
		log("Refreshing plotframe")
		self.plotlabels = []
		self.plotimages = []
		for index, plot in enumerate(self.plots):
			if plot == None:
				continue
			log("Adding plot")
			newplotimage = GUI.GUIHelper.cleanImage(plot, self.plotresx, self.plotresy)
			newplotlabel = Label(self, image=newplotimage)
			newplotlabel.grid(row=index+4, column=0)
			self.plotlabels.append(newplotlabel)
			self.plotimages.append(newplotimage)

	def addCarbonDensity(self):
		self.addPlot(Controller.get_carbon_density(currentMol(),traj=currentTraj()),0)

	def addTangentCorrelation(self):
		self.addPlot(Controller.get_tangent_correlation(currentMol(),traj=currentTraj()),0)

	def addPairCorrelation(self):
		self.addPlot(Controller.get_pair_correlation(currentMol(),traj=currentTraj()),0)

	def addPlot(self,path,index):
		log("Adding plot")
		self.plots[index]=path
		self.refresh()


	def update(self):
		log("plotframe updated")



statusframe = listBoxFrame(window, "Molecules")
statusframe.grid(column=0, row=0)
trajselectframe = listBoxFrame(window, "Traj Select")
trajselectframe.grid(column=1, row=0)
viewframe = listBoxFrame(window, "View")
viewframe.grid(column=2, row=0)
updateQueue = Button(window, text = "Update All", font = ("Arial Bold", 15), command = lambda: Controller.updateQueue())
updateQueue.grid(column=0, row=1)
plotsframe = plotFrame(window, width = 400, height = 600, bg = "gray")
plotsframe.grid(column=4, row=0, rowspan=500)
pic0, pic1, pic2, pic3 = Label(window),Label(window),Label(window),Label(window)
pics = [pic0, pic1, pic2, pic3]

# printframe = Frame(master = window, width = 400, height = 600)
# printframe.grid(column=1, row=0)
mols = []
printSelections = []


photos = [None, None, None, None]

# creating logger
logger = Logger.Logger("GUI")
def log(s):
	logger.log(s)

	# print queue.queue[queue.molToQueue[mol]]

# refreshes Status Frame
def refreshStatusFrame():
	statusframe.clear()
	for molclass, mol, milestone, status, line in Controller.queue.queue:
		statusframe.insert(END, "{!s} | {!s} {!s}".format(mol, milestone, status), mol)

# refreshes Traf Frame

def refreshTrajFrame():
	log("Refreshing Traj Frame")
	selectedMol = statusframe.selection()
	if not (selectedMol == None):
		options = [(selectedMol + " Traj"+str(i),i) for i in range(2,7)]
		options.insert(0,(selectedMol + " Init",1))
		log("options are: " + str(options))
		trajselectframe.refreshList(options)
		log("Selected molecule: " + str(statusframe.selection()))

def refreshViewFrame():
	log("Refreshing View Frame")
	options = [("Print" + str(i),("Print",i)) for i in range(1,5)]
	options.extend(("Video" + str(i),("Video",i)) for i in range(1,2))
	viewframe.refreshList(options)

def refreshPicture(index):
	pic = pics[index]
	pic.grid(column = 0, row=6, columnspan=3)

# changes the photo in photo[index]
def changePhoto(path, index=0):
	log("Changing Photo to: " + str(path))
	global photos
	maxwidth, maxheight = photos_resolution[index]
	photos[index] = GUI.GUIHelper.cleanImage(path, maxwidth, maxheight)
	pics[index].destroy()
	pics[index]=Label(window, image=photos[index])
	refreshPicture(index)

def playGif(path, index=0):
	log("Playing Gif")
	maxwidth, maxheight = photos_resolution[index]
	gif = GUI.GUIHelper.AnimatedGIF(window, path, maxwidth=maxwidth, maxheight=maxheight)
	pics[index].destroy()
	pics[index]=gif
	refreshPicture(index)

def currentMol():
	return statusframe.selection()

def currentMolClass():
	return Controller.queue.getClass(currentMol())

def currentTraj():
	return trajselectframe.selection()

def currentView():
	return viewframe.selection()

def openView():
	mol, traj = currentMol(), currentTraj()
	viewtype, viewnum = currentView()
	# open prints within the GUI
	log("Opening: " + str([mol, traj, viewtype, viewnum]))
	if (viewtype=="Print"):
		changePhoto(Controller.get_print(mol, traj=traj, view=viewnum))
	# open gifs otherwise
	else:
		playGif(Controller.get_video(currentMol(), traj=currentTraj(), view=viewnum))

statusframe.setOnselect(refreshTrajFrame)
trajselectframe.setOnselect(refreshViewFrame)
viewframe.setOnselect(openView)

refreshStatusFrame()
window.geometry(RESOLUTION)
photos_resolution = [(600.,400.),(250.,250.),(250.,250.),(250.,250.)]


window.mainloop()