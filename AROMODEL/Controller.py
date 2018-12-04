# import GUI
import Logger
import Queue
import Configure
import os
import ovitoHandler
import GUI.AnimatedGIF as gif
# Constructing Logger
logger = Logger.Logger("Controller")
queue = Queue.Queue()

def log(s):
	logger.log(s)


# updates entire queue
def updateQueue():
	for mol in queue.mols:
		updateMol(mol)

# updates a specific mol
def updateMol(mol):
	log("Updating mol: " + mol)
	queue.update(mol)

# shows a latest video given a molecules' milestone
def get_latest_video(molclass, mol, view=1):
	log("Showing latest video " + str([molclass, mol]))
	latesttraj = queue.get_latest_completed_trajectory(mol)
	if latesttraj==1:
		return log("No trajectories completed beyond initialization :(")
	return get_video(molclass, mol, latesttraj, view)

# shows a latest print given a molecules' milestone
def get_latest_print(molclass, mol, view=1):
	log("Showing latest print " + str([molclass, mol]))
	latesttraj = queue.get_latest_completed_trajectory(mol)
	return get_print(molclass, mol, latesttraj, view)

# def get_default_plot(mol, plot):
# 	return get_plot(queue.getClass(mol))

# shows latest carbon density
def get_carbon_density(mol, molclass=None, traj=None):
	if molclass is None:
		molclass = queue.getClass(mol)
	if traj is None:
		traj = queue.get_latest_completed_trajectory(mol)
	return get_plot(molclass, mol, traj, "CarbonDensity")

# shows latest tangent correlation
def get_tangent_correlation(mol, molclass=None, traj=None):
	if (molclass == None):
		molclass=queue.getClass(mol)
	if (traj==None):
		traj = queue.get_latest_completed_trajectory(mol)
	return get_plot(molclass, mol, traj, "TangentCorrelation")

# shows latest carbon density
def get_pair_correlation(mol, molclass=None, traj=None):
	if molclass is None:
		molclass = queue.getClass(mol)
	if traj is None:
		traj = queue.get_latest_completed_trajectory(mol)
	return get_plot(molclass, mol, traj, "PairCorrelation")

# general show functions
def get_plot(molclass, mol, traj, plot):
	log(("Showing plot ") + str([molclass, mol, traj, plot]))
	if (traj == 1):
		return
	traj = str(traj)
	trajplot = Configure.molGraphsPath(molclass,mol)+mol+"_"+plot+".png"
	log("Trajplot is: " + trajplot)
	if not os.path.exists(trajplot):
		log("Creating plot...")
		os.system("python " + Configure.Aromodel_Path +"/Analyze_SAM_Traj.py " +molclass + " " + mol + " " + traj)
	return trajplot

def get_video(mol, molclass=None, traj=None, view=1):
	if (molclass == None):
		molclass=queue.getClass(mol)
	if (traj == None):
		traj=queue.get_latest_completed_trajectory(mol)
	log("Showing video " + str([molclass, mol, traj, view]))
	if (traj == 1):
		traj = "INIT"
	traj = str(traj)
	trajvideo = Configure.molVideosPath(molclass, mol)+ mol + "_" + traj + "_view" + str(view) + ".gif"
	log("Trajvideo is: " + trajvideo)
	if not os.path.exists(trajvideo):
		inputfile = Configure.molResultsPath(molclass, mol)
		if (traj == "INIT"):
			inputfile = inputfile + "data." + mol
		else:
			inputfile = inputfile + mol + "_" + traj + ".lammpstrj"
		ovitoHandler.videoView(inputfile, view)
	return trajvideo

def get_print(mol, molclass=None, traj=None, view=1):
	
	if (molclass == None):
		molclass=queue.getClass(mol)
	if (traj == None):
		traj=queue.get_latest_completed_trajectory(mol)
	log("Showing print " + str([molclass, mol, traj, view]))
	traj = str(traj)
	if (traj == "Last"):
		return get_latest_print(molclass, mol, view)
	if (traj == "1"):
		traj = "INIT"
	traj = str(traj)
	trajprint = Configure.molPrintsPath(molclass, mol)+ mol + "_" + traj + "_view" + str(view) + ".png"
	if not os.path.exists(trajprint):
		inputfile = Configure.molResultsPath(molclass, mol)
		if (traj == "INIT"):
			inputfile = inputfile + "data." + mol
		else:
			inputfile = inputfile + mol + "_" + traj + ".lammpstrj"
		ovitoHandler.printView(inputfile, view)
	return trajprint

if __name__=='__main__': main()
