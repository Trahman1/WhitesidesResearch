import Configure
import Logger
import os

# class containing current Queue
# queuelines is a line by line of the Queue.txt
# queue is an array containing [molclass, mol, milestone, status, linenumber]
# molToQueue is a mapping from the mol name to its queue entry
logger = Logger.Logger("Queue")

def log(s):
	logger.log(s)

class Queue():
	def __init__(self):
		self.queuelines = []
		self.queue = []
		self.molToQueue = {}
		self.mols = []
		self.initQueue()
		# print self.molToQueue

	def initQueue(self):
		log("Opening Queue with queuepath: " + Configure.Queue_Path)
		with open(Configure.Queue_Path, 'r') as queuefile:
			self.queuelines = queuefile.readlines()
			self.queue = self.queuelines[:]
		
		# split elements in a row as separate elements in a list
		# row contains molclass, mol, 
		for index, line in enumerate(self.queue):
				self.queue[index] = line.strip("\n").split()+[index]
				# print self.queue[index]
		
		# remove blank lines from queue
		self.queue = filter(lambda line: len(line)>1, self.queue)
		
		# init molTOQueue and mols
		for index, line in enumerate(self.queue):
			self.molToQueue[(self.queue[index][1])]=index
			self.mols.append(self.queue[index][1])
		# print self.queue

	def getRow(self, mol):
		# molclass, mol, milestone, status, line
		return self.queue[self.molToQueue[mol]]

	def getClass(self, mol):
		molclass, mol, milestone, status, line = self.getRow(mol)
		return molclass

	def writeQueue(self):
		# print queuelines
		newfile = open(Configure.Queue_Path, 'w')
		newfile.writelines(self.queuelines)
		newfile.close()
		log("Refreshed Queue")

	def printQueue(self):
		print self.queue
		for molclass, mol, milestone, status, linenumber in self.queue:
			print ("{!s}:{!s} is {!s} {!s}").format(molclass, mol, milestone, status)

	def update(self, mol):
		molclass, mol, milestone, status, line = self.queue[self.molToQueue[mol]]
		milestone, status = self.getStatus(mol)
		self.queue[self.molToQueue[mol]] = molclass, mol, milestone, status, line
		self.queuelines[line] = molclass + " " + mol + " " + milestone + " " + status + "\n"
		log("Status of " + molclass + ":" + mol +" is " + milestone + "_" + status)
		self.writeQueue()

	def updateAll(self):
		log("Updating All")
		for mol in self.mols:
			self.update(mol)

	# function that returns latest completed trajectory
	def get_latest_completed_trajectory(self, mol):
		molclass, mol, milestone, status, line = self.getRow(mol)
		if (milestone == "INIT"):
			return 1
		latesttraj = int(milestone[-1:])
		if ((status == "RUNNING") or (status == "AWAITING")):
			latesttraj = latesttraj-1
		return latesttraj
		# print "Queue read!"
		# print "Queue updated!"
		# print self.queuelines
		# print queue
		# for index, line in enumerate(queue):
		# 	queue[index] = [index] + queue[index] 


	# returns the status in form of (TRAJX, RUNNING)
	def getStatus(self, mol):
		log("Getting status of: " + mol + "...")
		statusFiles = self.remainingStatusFiles(mol)
		statusstr = None
		for milestone, startmarker, finishmarker in statusFiles:
			if (checkFile(mol,finishmarker)):
				continue
			elif(checkFile(mol,startmarker)):
				statusstr = (milestone,"RUNNING")
				break
			else:
				statusstr = (milestone,"AWAITING")
				break
		if (statusstr == None):
			statusstr = ("TRAJ6", "COMPLETE")
		log("Status for " + mol + ": " + statusstr[0]+"_"+statusstr[1])
		return statusstr

	# helper function for getStatus - creates list of files to check
	def remainingStatusFiles(self, mol):
		molclass, mol, milestone, status, line = self.getRow(mol)

		latestMilestone = 1
		if (milestone[:4]=="TRAJ"):
			latestMilestone = int(milestone[4:])
		milestones = ["TRAJ" + str(i) for i in range(2,7)]
		statusFiles = []

		if (latestMilestone == 1):



			statusFiles.append(("INIT", "Init_"+mol+".lammpstrj", "restart."+mol+"_init_1"))
			latestMilestone = latestMilestone+1
		for i in range(latestMilestone, 7):
			temp = "400"
			if i >= 4:
				temp = "300"
			statusFiles.append(("TRAJ"+str(i), mol+"_"+str(i)+".lammpstrj", "restart."+mol+"_"+str(temp)+"_"+str(i)))
		return statusFiles
	



	# status, startmarker, finishmarker

def checkFile(mol, file):
	log("Checking for " + mol + ": " + file)
	filepath = Configure.Comet_Path.strip("%s")+mol+"/"+file
	call = "ssh " + Configure.Comet_Login + " [ -f " + filepath + " ] && echo 'exist' || echo 'notExist'"
	status = os.popen(call).read()
	if status.strip("\n")=='exist':
		log("Found " + mol + ": " + file)
	else:
		log("Did not find " + mol + ": " + file)
	return (status.strip("\n")=='exist')

def main():
	# updateQueue()
	queue = Queue()
	queue.updateAll()

if __name__=='__main__': main()

# status("OEG1")
# status("OEG2")
# status("OEG3")
# status("OEG4")
# status("OEG5")
# status("OEG6")
# status("OEG8")




