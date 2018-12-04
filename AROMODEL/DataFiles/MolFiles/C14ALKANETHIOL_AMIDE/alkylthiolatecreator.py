# function that creates OEGx in current directory using openBabel
import sys
import os

groupname = "C14ALKANETHIOL_AMIDE"

def createMol(numGroup, rectifyingGroup):
	bashcall = "obabel -:\"[S]"
	repeat = "C"
	chemname = groupname
	
	# if numGroup = -1, creates control group
	if (numGroup==-1):
		chemname += "_control"
	else:
		chemname += str(numGroup)
	filename = chemname+".xyz"
	backbonesize = 12
	for index in range(backbonesize):
		if (index == numGroup):
			bashcall+= rectifyingGroup
		else:
			bashcall+=repeat
	
	# if numGroup = -1, creates control group
	if (numGroup==(-1)):
		bashcall += repeat
	bashcall+="C\" -O " + filename + " --gen3d"
	print "bashcall is: " + bashcall
	os.system(bashcall)
	return chemname

def main():
	rectifyingGroup = "C(=O)N"
	if (len(sys.argv)==1):
		print "No argument given; Constructing groups 1 to 10 for : " + groupname
		for n in range(9):
			dirname = createMol(n+1, rectifyingGroup)
			# os.system("mkdir " + dirname)
			# os.system("cp Aromodel.py " + dirname)
			# os.system("cp gold.data " + dirname )
			# os.system("mv " + dirname + ".xyz " + dirname)
	else:
		numO = int(sys.argv[1])
		print "Constructing " + groupname + str(numO)
		createMol(numO, rectifyingGroup)

if __name__=='__main__': main()

