# function that creates OEGx in current directory using openBabel
import sys
import os
import numpy as np

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
	rotateMolXtoZ(filename)
	return chemname

# constructs rotation matrix that rotates prevector to the postvector
def rotationMatrix(prevector, postvector):
	# cross = np.cross(vector1, vector2)
	# x,y,z = cross
	print prevector, postvector
	prevector, postvector = list(map(lambda v: v/np.linalg.norm(v), [prevector,postvector]))
	print prevector, postvector
	c = np.dot(prevector, postvector)
	s = np.cross(prevector, postvector)
	crossmatrix = np.array([[0,-s[2],s[1]],[s[2],0,-s[0]],[-s[1],s[0],0]])
	rotationmatrix = np.identity(3) + crossmatrix + 1./(1.+c)*np.dot(crossmatrix,crossmatrix)
	return rotationmatrix


def writeXYZline(atom, pos):
	x, y, z = pos[0], pos[1], pos[2]
	return "{atom} {x} {y} {z} \n".format(atom = atom, x=x, y=y, z=z)

def rotateMolXtoZ(filename):
	with open(filename) as f:
		lines = f.readlines()
	Spos = np.asarray(list(map(lambda x: float(x), lines[2].split()[1:])))
	Cpos = np.asarray(list(map(lambda x: float(x), lines[3].split()[1:])))
	backbone = Cpos-Spos

	rotation = rotationMatrix(backbone, np.asarray([0,1,0]))
	# print lines
	with open(filename, "w") as f:
		f.write(lines[0])
		f.write(lines[1])
		for line in lines[2:]:
			line = line.split()
			atom = line[0]
			pos = np.asarray(list(map(lambda x: float(x), line[1:])))
			pos = np.dot(rotation, pos)
			# print writeXYZline(atom, pos)
			f.write(writeXYZline(atom, pos))
	
	


def main():
	if (len(sys.argv)==1):
		print "No argument given; Constructing OEG1 - OEG 10"
		for n in range(9):
			chemname = createMol(n+1)
			# os.system("mkdir " + dirname)
			# os.system("cp Aromodel.py " + dirname)
			# os.system("cp gold.data " + dirname )
			# os.system("mv " + dirname + ".xyz " + dirname)
	else:
		numO = int(sys.argv[1])
		print "Constructing OEG" + str(numO)
		filename = createMol(numO) +".xyz"

if __name__=='__main__': main()

