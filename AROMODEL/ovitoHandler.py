# Handler python module (written in Python 2) for ovito files

import Configure
import os
import Logger
import sys

logger = Logger.Logger("OvitoHandler")

def log(s):
	logger.log(s)

def createPath(path):
	if not os.path.exists(path):
			os.mkdir(path)



def printView(inputfile, view, outputfile=None, renderer = "t"):
	log("Printing view: " + str([inputfile, view, outputfile, renderer]))
	if (outputfile==None):
		printpath = (os.path.dirname(inputfile))+"/Prints/"
		createPath(printpath)
		outputname = os.path.basename(inputfile).split(".")
		if outputname[0]=="data":
			outputname = outputname[1]+"_INIT"
		else:
			outputname=outputname[0]
		outputfile = printpath + outputname + "_view" + str(view) + ".png"
	viewfile = Configure.Ovito_Views_Path + "View" + str(view) + ".ovito"
	printscript = Configure.Ovito_Scripts_Path + "print.py" 
	bashcall = Configure.Ovitos_Path + " -o " + viewfile + " " + printscript + " " + inputfile + " printtemplated " + renderer + " p " + outputfile
	log(bashcall)
	os.system(bashcall)

def videoView(inputfile, view, outputfile=None, renderer = "t"):
	log("Video-ing view: " + str([inputfile, view, outputfile, renderer]))
	if (outputfile==None):
		videopath = (os.path.dirname(inputfile))+"/Videos/"
		createPath(videopath)
		outputname = os.path.basename(inputfile).split(".")
		if outputname[0]=="data":
			outputname = outputname[1]+"_INIT"
		else:
			outputname=outputname[0]
		outputfile = videopath + outputname + "_view" + str(view) + ".gif"
	viewfile = "video1"
	script = Configure.Ovito_Scripts_Path + "print.py" 
	bashcall = Configure.Ovitos_Path + " " + script + " " + inputfile + " " + viewfile + " " + renderer + " v " + outputfile
	log(bashcall)
	os.system(bashcall)

def main():
	script, inputfile = sys.argv
	for view in [1,2,3,4]:
		printView(inputfile, view, renderer = "t")


if __name__=='__main__': main()

