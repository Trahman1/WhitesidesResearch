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



# Remote Paths
Comet_Login = "trahman2@comet.sdsc.edu"
Comet_Path = "/oasis/scratch/comet/trahman2/temp_project/Aromodel/%s" # % Directory_Name
Orca_Path = "/oasis/scratch/comet/cjpais/temp_project/programs/orca_3_0_3_linux_x86-64/orca"



c2c = "scp %s " + Comet_Login + ":" + Comet_Path # % (File, Directory_Name)
c2l = "scp " + Comet_Login + ":" + Comet_Path + "/%s ./" # % (Directory_Name, File)
SBATCH = "sbatch " + Comet_Path + "/%s" # %( Directory_Name, Submit_Script)

#cd /oasis/scratch/comet/trahman2/temp_project/