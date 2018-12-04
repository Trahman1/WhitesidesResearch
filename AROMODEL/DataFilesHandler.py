# !/usr/bin/python

import Configure 

# returns surface file and returns .data
def getSurface(surfaceName):
	if (surfaceName == "gold"):
		return Configure.SurfaceFiles_Path + "gold.data"
	else:
		return Configure.SurfaceFiles_Path + surfaceName + ".data" 

def getGoldSurface():
	return getSurface("gold")

# returns mol file and appends .xyz
def getMol(molName, molClass):
	return Configure.MolFiles_Path + molClass + "/" + molName + ".xyz"

def getTraj(molName,molClass,TrajNumber):
	return Configure.Results_Path + molClass + "/" + molName + "/" + molName + "_" + TrajNumber + ".lammpstrj"