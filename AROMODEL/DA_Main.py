import Atom
import Molecule
import OPLS
import System
import sys
import DA_Polymer
import Surface
import glob
import numpy as np
import os

# This script is specific to the DA project, it takes in lammps data files
# and generates a polymer system

def main():
    Script, Mult = sys.argv
    Density = 0.01
    Mult = int(Mult)
    File_List = glob.glob('*.data')
    print len(File_List)
    print File_List
    Mol_Temp_List = []
    Mol_Temp_List.append(Surface.Surface('alkane-monolayer_2.lammps'))
    i = 0
    Total_Mass = 0.0
    MW = 0.0
    for File in File_List:
        Mol_Temp_List.append(DA_Polymer.DA_Polymer(File))
        Total_Mass += Mol_Temp_List[i].MW
        MW = Mol_Temp_List[i].MW
        i += 1
    print i
    Avogadro = 6.0221413e23
    Comp_List = np.ones(i+2, dtype=int)
    Comp_List = Comp_List*Mult
    Comp_List[0] = 1
    Total_Mols = float(i*Mult + 1)/Avogadro
    Molar_Volume = MW/Density
    Volume = Molar_Volume*Total_Mols

    #Name = File_List[0].split('.')[1].split('_')[0] + "_%s" % i*Mult
    
    Name = os.getcwd().split('/')[-1] + "_%s" % int(i*Mult)
    print Name
    DA_System = System.System(Mol_Temp_List, Comp_List, Box_Length, Name)
    DA_System.Gen_Rand()
    DA_System.Write_LAMMPS_Data()
    DA_System.Run_Lammps_Init()
    
    DA_System.Temperature = 800
    DA_System.Run_Lammps_NPT(GPU=True)
    DA_System.Run_Lammps_NPT(GPU=True)
    #System.Run_Glass_Transition(DA_System, -100, Ramp_Steps = 1000000, Equil_Steps = 1000000, T_End = 600)
    DA_System.Run_Lammps_NPT(GPU=True)
    DA_System.Run_Lammps_NPT(GPU=True)
    DA_System.Run_Lammps_NPT(GPU=True)
    DA_System.Temperature = 600
    DA_System.Run_Lammps_NPT(GPU=True)
    DA_System.Run_Lammps_NPT(GPU=True)
    DA_System.Run_Lammps_NPT(GPU=True)
    DA_System.Run_Lammps_NPT(GPU=True)
    DA_System.Run_Lammps_NPT(GPU=True)
    DA_System.Run_Lammps_NPT(GPU=True)
    System.Run_Glass_Transition(DA_System, 20, Ramp_Steps = 100000, Equil_Steps = 100000, T_End = 100)


if __name__=='__main__': main()


