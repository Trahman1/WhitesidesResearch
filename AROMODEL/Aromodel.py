#!/usr/bin/python

import Atom
from Configure import Results_Path as resPath
from Configure import Aromodel_Path as aroPath
import Molecule
import OPLS
import System
import sys
import Surface
import SAM
import time
import glob
import numpy as np
import os
import glob
import DataFilesHandler as data

# molPath = ""
# surfacePath = ""

# def read_input(dirname):
#     f = open(glob.glob(dirname+'/*.aro.in')[0],"r")
#     molFile=f.readline()
#     surfaceFile = f.readline()
#     Name = f.readline()
#     return [molFile, surfaceFile, Name]

def main():
    # print sys.argv
    try:
        scriptName, surface, molClass, mol= sys.argv
        cmd = "ttab python " + aroPath + "/Aromodel.py " + surface + " " + molClass + " " + mol + " RunningTab"
        print cmd
        os.system(cmd)
    except:
        scriptName, surface, molClass, mol, thistab = sys.argv
        runAll(mol, surface, molClass)
    # try:
    #     print "trying"
    #     scriptName, surface, molClass, mol= sys.argv
    #     print "tried!"
    #     cmd = "ttab python " + Configure.Aromodel_Path + "/Aromodel.py " + surface + " " + molClass + " " + mol + " RunningTab"
    #     print cmd
    #     os.system(cmd)
    # except:
    #     scriptName, surface, molClass, mol, thistab = sys.argv
    #     run(mol, surface, molClass)
    # molFile = 'OEG2.xyz'
    # surfaceFile = 'gold.data'
    # Name = 'OEG2'
    # print (molFile, surfaceFile, Name)

def resultPath(mol, molClass):
    target = resPath + molClass + "/" + mol + "/"
    if (not (os.path.exists(target))):
        print "making dir: " + target   
        os.mkdir(target)
    # else:
    print "changing dir: " + target
    os.chdir(target)
    # return target


def runAll(mol, surface, molClass):
    resultPath(mol,molClass)
    au_a = 4.06 # Lattice constant of gold
    relative_coverage = 1.0
    Density = 21.6*relative_coverage # Angstroms^2/chain
    GOLD_Obj = Surface.Surface(data.getSurface(surface)) # Make gold "Surface Object"
    Mol_Obj = Molecule.Molecule(mol, data.getMol(mol, molClass)) # Make Molecular objects
    Mol_Obj.UnConverged = True
    Box_Length = GOLD_Obj.Box_Length
    Area = Box_Length[0]*Box_Length[1]
    Num_SAMs = int(np.sqrt(Area/Density))
    Mol_Obj.Set_Up_FF(run_orca=True, local = True)# Parameterize Molecule Object
    OPLS.Assign_OPLS(Mol_Obj, ChelpG = False) # Parameterize Molecule Object
    OPLS.Assign_OPLS(GOLD_Obj, ChelpG = False) # Parametrize Au slab
    SAM_System = SAM.SAM(Mol_Obj, GOLD_Obj, Num_SAMs, Box_Length, mol)
    SAM_System.Gen_SAM()
    SAM_System.Write_LAMMPS_Data()
    SAM_System.Run_Lammps_Init()
    SAM_System.Run_Lammps_NPT(Temp_Out = 400, time_steps = 5000000)
    SAM_System.Run_Lammps_NPT(Temp_Out = 400, time_steps = 5000000)
    SAM_System.Run_Lammps_NPT(Temp_Out = 300, time_steps = 5000000)
    SAM_System.Run_Lammps_NPT(Temp_Out = 300, time_steps = 5000000)
    SAM_System.Run_Lammps_NPT(Temp_Out = 300, time_steps = 5000000)
    print "Aromodel.py complete for: " + mol + " molecule on " + surface + " surface!"

    #BHJ.Run_Lammps_Init(Nodes=2)
    #BHJ.Run_Lammps_NPT(Nodes=2)
    #BHJ.Run_Lammps_NPT(Nodes=2)
    #System.Run_Glass_Transition(BHJ, 25, Ramp_Steps = 100000, Equil_Steps = 100000, T_End = 100, Nodes=1)

if __name__=='__main__': main()


