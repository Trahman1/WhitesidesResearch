# function that creates OEGx in current directory using openBabel
import sys
import os

groupname = "C12ALKANETHIOL"

def createMol(numGroup, rectifyingGroup):
	bashcall = "obabel -:\"[S]"
	repeat = "C"
	if (rectifyingGroup=="control"):
		filename = "C12ALKANETHIOL"
		bashcall+="CCCCCCCCCCCC\" -O " + filename + " --gen3d"

	elif (rectifyingGroup=="thiourea"):
		rectifyingGroup = "NC(=S)N"
		rectifierName = "NHCSNH"
		backbonesize = 10
		chemname = "C" + str(int(numGroup)-1) + "-" + rectifierName + "-" + "C" + str(10-int(numGroup))
		filename = chemname + ".xyz"
		for index in range(1,backbonesize):
			if (index == int(numGroup)):
				print "rectifyingGroup is " + rectifyingGroup
				bashcall+= rectifyingGroup
			else:
				bashcall+=repeat
		bashcall+="C\" -O " + filename + " --gen3d"
		print "bashcall is: " + bashcall
		os.system(bashcall)

	elif (rectifyingGroup=="urea"):
		rectifyingGroup = "NC(=O)N"
		rectifierName = "NHCONH"
		backbonesize = 10
		chemname = "C" + str(int(numGroup)-1) + "-" + rectifierName + "-" + "C" + str(10-int(numGroup))
		filename = chemname + ".xyz"
		for index in range(1,backbonesize):
			if (index == int(numGroup)):
				print "rectifyingGroup is " + rectifyingGroup
				bashcall+= rectifyingGroup
			else:
				bashcall+=repeat
		bashcall+="C\" -O " + filename + " --gen3d"
		print "bashcall is: " + bashcall
		os.system(bashcall)

	elif (rectifyingGroup=="amide"):
		rectifyingGroup = "CC(=O)N"
		rectifierName = "CONH"
		backbonesize = 10
		chemname = "C" + str(numGroup) + "-" + rectifierName + "-" + "C" + str(10-int(numGroup))
		filename = chemname + ".xyz"
		for index in range(1,backbonesize):
			if (index == int(numGroup)):
				print "rectifyingGroup is " + rectifyingGroup
				bashcall+= rectifyingGroup
			else:
				bashcall+=repeat
		bashcall+="C\" -O " + filename + " --gen3d"
		print "bashcall is: " + bashcall
		os.system(bashcall)
		print "printing second amide"
		rectifyingGroup = "NC(=O)C"
		rectifierName = "NHCO"
		bashcall = "obabel -:\"[S]"
		print rectifierName
		chemname = "C" + str(int(numGroup)-1) + "-" + rectifierName + "-" + "C" + str(11-int(numGroup))
		filename = chemname + ".xyz"
		backbonesize = 10
		for index in range(1,backbonesize):
			if (index == int(numGroup)):
				print "rectifyingGroup is " + rectifyingGroup
				bashcall+= rectifyingGroup
			else:
				bashcall+=repeat
		bashcall+="C\" -O " + filename + " --gen3d"
		print "bashcall is: " + bashcall
		os.system(bashcall)

def main():
	# rectifyingGroup = "NC(=O)C"
	# rectifierName = "NHCO"
	if (len(sys.argv)==3):
		name, rectifyingGroup, numGroup = sys.argv
		createMol(numGroup, rectifyingGroup)
	else:
		name, numGroup = sys.argv
		for rectifyingGroup in ["urea", "thiourea", "amide"]:
			createMol(numGroup, rectifyingGroup)

	# else:
	# 	name, numGroup = sys.argv
	# 	createMol(numGroup, "CC(=O)N", "CONH")
	# 	createMol(numGroup, "NC(=O)C", "NHCO")
	# 	createMol(numGroup, "NC(=O)N", "NHCONH")
	# 	createMol(numGroup, "NC(=S)N", "NHCSNH")



if __name__=='__main__': main()

