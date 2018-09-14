# !/usr/bin/python

from Configure import MolFiles_Path as mfPath
from Configure import SurfaceFiles_Path as sfPath

# returns surface file and returns .data
def getSurface(surfaceName):
	if (surfaceName == "gold"):
		return sfPath + "gold.data"
	else:
		return sfPath + surfaceName + ".data" 

def getGoldSurface():
	return getSurface("gold")

# returns mol file and appends .xyz
def getMol(molName):
	if (molName[:3] == "OEG"):
		return mfPath + "OEG/" + molName + ".xyz"
	else:
		return mfPath + molName + ".xyz"
