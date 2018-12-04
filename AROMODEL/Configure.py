#! usr/bin/python

"""
    This file contains paths to different executables 
    and template files, so that this software can be 
    configured easily on different machines.
"""


# Local Paths
Research_Path = "/Users/tamjidrahman/Desktop/Whitesides_Research"
Aromodel_Path = Research_Path + "/AROMODEL"
Template_Path = Aromodel_Path + "/Templates/"
DataFiles_Path = Aromodel_Path + "/DataFiles/"
SurfaceFiles_Path = DataFiles_Path + "SurfaceFiles/"
MolFiles_Path = DataFiles_Path + "MolFiles/"
Results_Path = Research_Path + "/Results/"
Queue_Path = Aromodel_Path + "/Queue"
Ovito_Path = Aromodel_Path + "/Ovito/"
Ovito_Views_Path = Ovito_Path+"Views/"
Ovito_Scripts_Path = Ovito_Path+"Scripts/"
Ovito_Print_Path = Ovito_Scripts_Path+"print.py"
Ovitos_Path = "~/../../Applications/Ovito.app/Contents/MacOS/ovitos"

def molResultsPath(mol_class, mol):
	return Results_Path + mol_class + "/" + mol + "/"

def molGraphsPath(mol_class, mol):
	return molResultsPath(mol_class, mol) + "Graphs/"

def molPrintsPath(mol_class, mol):
	return molResultsPath(mol_class, mol) + "Prints/"

def molVideosPath(mol_class, mol):
	return molResultsPath(mol_class, mol) + "Videos/"


# Remote Paths
Comet_Login = "trahman2@comet.sdsc.edu"
Comet_Path = "/oasis/scratch/comet/trahman2/temp_project/Aromodel/%s" # % Directory_Name
Orca_Path = "/oasis/scratch/comet/cjpais/temp_project/programs/orca_3_0_3_linux_x86-64/orca"



c2c = "scp %s " + Comet_Login + ":" + Comet_Path # % (File, Directory_Name)
c2l = "scp " + Comet_Login + ":" + Comet_Path + "/%s ./" # % (Directory_Name, File)
SBATCH = "sbatch " + Comet_Path + "/%s" # %( Directory_Name, Submit_Script)

#cd /oasis/scratch/comet/trahman2/temp_project/