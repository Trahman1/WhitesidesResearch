# function that creates OEGx in current directory using openBabel
import sys
import os

def createMol(numO):
	bashcall = "obabel -:\"[S]"
	repeat = "CCO"
	terminate = "C"
	chemname = "OEG" + str(numO)
	filename = chemname+".xyz"
	for ocount in range(numO):
		bashcall+=repeat
	bashcall+="C\" -O " + filename + " --gen3d"
	os.system(bashcall)
	return chemname

def main():
	if (len(sys.argv)==1):
		print "No argument given; Constructing OEG1 - OEG 10"
		for n in range(9):
			dirname = createMol(n+1)
			# os.system("mkdir " + dirname)
			# os.system("cp Aromodel.py " + dirname)
			# os.system("cp gold.data " + dirname )
			# os.system("mv " + dirname + ".xyz " + dirname)
	else:
		numO = int(sys.argv[1])
		print "Constructing OEG" + str(numO)
		createMol(numO)

if __name__=='__main__': main()

