#!/usr/bin/python

import Atom
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

# molPath = ""
# surfacePath = ""

# def read_input(dirname):
#     f = open(glob.glob(dirname+'/*.aro.in')[0],"r")
#     molFile=f.readline()
#     surfaceFile = f.readline()
#     Name = f.readline()
#     return [molFile, surfaceFile, Name]

def main():
    molFile = 'OEG8.xyz'
    surfaceFile = 'gold.data'
    Name = 'OEG8'
    run(molFile, surfaceFile, Name)


def run(molFile, surfaceFile, Name):
    au_a = 4.06 # Lattice constant of gold
    relative_coverage = 1.0
    Density = 21.6*relative_coverage # Angstroms^2/chain
    GOLD_Obj = Surface.Surface(surfaceFile) # Make gold "Surface Object"
    SAM_Obj = Molecule.Molecule(molFile) # Make Molecular objects
    SAM_Obj.UnConverged = True
    Box_Length = GOLD_Obj.Box_Length
    Area = Box_Length[0]*Box_Length[1]
    Num_SAMs = int(np.sqrt(Area/Density))
    SAM_Obj.Set_Up_FF(run_orca=True, local = True)# Parameterize Molecule Object
    OPLS.Assign_OPLS(SAM_Obj, ChelpG = False) # Parameterize Molecule Object
    OPLS.Assign_OPLS(GOLD_Obj, ChelpG = False) # Parametrize Au slab
    SAM_System = SAM.SAM(SAM_Obj, GOLD_Obj, Num_SAMs, Box_Length, Name)
    SAM_System.Gen_SAM()
    SAM_System.Write_LAMMPS_Data()
    SAM_System.Run_Lammps_Init()
    SAM_System.Run_Lammps_NPT(Temp_Out = 400, time_steps = 5000000)
    SAM_System.Run_Lammps_NPT(Temp_Out = 400, time_steps = 5000000)
    SAM_System.Run_Lammps_NPT(Temp_Out = 300, time_steps = 5000000)
    SAM_System.Run_Lammps_NPT(Temp_Out = 300, time_steps = 5000000)
    SAM_System.Run_Lammps_NPT(Temp_Out = 300, time_steps = 5000000)

    #BHJ.Run_Lammps_Init(Nodes=2)
    #BHJ.Run_Lammps_NPT(Nodes=2)
    #BHJ.Run_Lammps_NPT(Nodes=2)
    #System.Run_Glass_Transition(BHJ, 25, Ramp_Steps = 100000, Equil_Steps = 100000, T_End = 100, Nodes=1)

if __name__=='__main__': main()


