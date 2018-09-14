#! usr/bin/python
import Molecule
import random
import numpy as np
from copy import deepcopy
import os
import subprocess
import time
import Configure
import Parallel
import pickle
import Bond
import dataReader
# Class defining an MD system for simulation of SAMS with LAMMPS

class SAM(object):
    """
        Class defining an MD system for simulation of SAMs with LAMMPS
        instance variables:
        Molecule_List = List of Molecule objects
        Composition_List = List of integers defining the number of each molecule type in the system
        Box_Size = float (Angstrom)
        Atom_Params = [[Mass, Sigma, Epsilon]] List of floats
        Bond_Params = [[KB, RO]] List of floats
        Angle_Params = [[KA, A0]] List of floats
        Dihedral_Params = [[V1,V2,V3,V4]] List of floats
        Improper_Params = [[KI, AI]] List of floats
        Num_Atoms = float
        
        """
    def __init__(self, SAM, Surface, Num_SAMs,Box_Length, Name, Solvent=False):
        self.Name = Name
        self.SAM = SAM
        self.Surface = Surface
        self.Num_SAMs = Num_SAMs
        self.Box_Length = Box_Length
        self.Molecule_List = []
        self.Current_Restart = ""
        self.Temperature = 400
        if Solvent != False:
            self.Solvent = Solvent
            self.Solvent_Log = True
            self.Atom_Params, self.Bond_Params, self.Angle_Params, self.Dihedral_Params, self.Improper_Params = Molecule.Assign_Lammps([SAM, Surface, Solvent])
        else:
            self.Solvent_Log = False
            self.Atom_Params, self.Bond_Params, self.Angle_Params, self.Dihedral_Params, self.Improper_Params = Molecule.Assign_Lammps([SAM, Surface])

        return

    def Gen_SAM(self):
        self.Molecule_List.append(self.Surface)
        N=0
        M=0
        j = 0
        for atom in self.SAM.Atom_List:
            j +=1
            if atom.Element == 'S':
                pos_s = atom.Position
                temp_atom = deepcopy(atom)
        i = 1
        for Molecule_Obj in self.Molecule_List:
            for Atom_Obj in Molecule_Obj.Atom_List:
                Atom_Obj.System_ID = i
                i += 1

        surface_list = []
        for atom in self.Surface.Atom_List:
            if atom.Position[2] >= 39.6:
                surface_list.append(atom.Position)
        Position = np.array([2.34,1.4, 9.0 + self.SAM.z_radius], dtype = float)
        for i in range(self.Num_SAMs):
            for j in range(self.Num_SAMs):
                M +=1
                temp_mol = deepcopy(self.SAM)
                temp_mol.Mol_ID = M
                temp_mol.COM += Position
                temp_mol.Adjust_COM(Random=False)
                self.Molecule_List.append(temp_mol)
                Position[0] += 4.06*(20./float(self.Num_SAMs))
            Position[0] = 2.34 + ((i+1)%2)*2.03*(20./float(self.Num_SAMs))
            Position[1] += 7.0321*(20./float(self.Num_SAMs))/2.0
        Num_Molecules = len(self.Molecule_List)
        
        if self.Solvent_Log:
            J = 1000
            for j in range(J):
                Temp_Mol = deepcopy(self.Solvent)
                Deposited = False
                while not Deposited:
                    Temp_Mol.COM = np.asarray([5 +random.random()*(self.Box_Length[0] - 15),5+ random.random()*(self.Box_Length[1] - 15), 30 + random.random()*(self.Box_Length[2]-50)], dtype = float)
                    M += 1
                    Temp_Mol.Mol_ID = M
                    Temp_Mol.Adjust_COM()
                    self.Molecule_List.append(Temp_Mol)
                    print "Depositing", j
                    Deposited = True

        return



    def Write_LAMMPS_Data(self, Dihedral = False):
        """
        Function for writing LAMMPS data file
        """
        if Dihedral:
            self.Data_File = "Dihedral.data"
        else:
            self.Data_File = "data." + self.Name
            File = open(self.Data_File, 'w')
            File.write('LAMMPS data file via SAM.Write_LAMMPS_Data()\n\n')

            # Find number of atoms, bonds, dihedrals, impropers in the system
            self.Num_Atoms = 0
            self.Num_Bonds = 0
            self.Num_Angles = 0
            self.Num_Dihedrals = 0
            self.Num_Impropers = 0
            i = 0
            for Mol_Obj in self.Molecule_List:
                self.Num_Atoms += len(Mol_Obj.Atom_List)
                self.Num_Bonds += len(Mol_Obj.Bond_List)
                self.Num_Angles += len(Mol_Obj.Angle_List)
                self.Num_Dihedrals += len(Mol_Obj.Dihedral_List)
                self.Num_Impropers += len(Mol_Obj.Improper_List)
                i += 1


            File.write('%d atoms\n' % self.Num_Atoms)
            File.write('%d atom types\n' % len(self.Atom_Params))
            File.write('%d bonds\n' % self.Num_Bonds)
            File.write('%d bond types\n' % len(self.Bond_Params))
            File.write('%d angles\n' % self.Num_Angles)
            File.write('%d angle types\n' % len(self.Angle_Params))
            if self.Num_Dihedrals > 0:
                File.write('%d dihedrals\n' % self.Num_Dihedrals)
                File.write('%d dihedral types\n' % len(self.Dihedral_Params))

            if self.Num_Impropers > 0:
                File.write('%d impropers\n' % self.Num_Impropers)
                File.write('%d improper types\n' % len(self.Improper_Params))

            try:
                if len(self.Box_Length) > 1:
                    File.write('\n\n0.0000 %.4f xlo xhi\n' % self.Box_Length[0])
                    File.write('0.0000 %.4f ylo yhi\n' % self.Box_Length[1])
                    File.write('0.0000 %.4f zlo zhi\n' % self.Box_Length[2])

            except:
                File.write('\n\n0.0000 %.4f xlo xhi\n' % self.Box_Length)
                File.write('0.0000 %.4f ylo yhi\n' % self.Box_Length)
                File.write('0.0000 %.4f zlo zhi\n' % self.Box_Length)

            File.write('\n\nMasses\n\n')
            i = 1
            for Params in self.Atom_Params:
                File.write('%d %.3f\n' % ( i, Params[0]))
                i += 1

            File.write('\n\nPair Coeffs # lj/cut/coul/long\n\n')
            i = 1
            for Params in self.Atom_Params:
                File.write('%d %.3f %.3f\n' % (i, Params[2], Params[1]))
                i += 1

            File.write('\n\nBond Coeffs # harmonic\n\n')
            i = 1
            for Params in self.Bond_Params:
                File.write('%d %.4f %.4f\n' % (i, Params[0], Params[1])) # Its possible that im missing a factor of 2 here
                i += 1

            File.write('\n\nAngle Coeffs # harmonic\n\n')
            i = 1
            for Params in self.Angle_Params:
                File.write('%d %.4f %.4f\n' % (i, Params[0], Params[1])) # Its possible that im missing a factor of 2 here
                i += 1
            if self.Num_Dihedrals > 0:
                File.write('\n\nDihedral Coeffs # opls\n\n')
                i = 1
                for Params in self.Dihedral_Params:
                    File.write('%d %.4f %.4f %.4f %.4f\n' % (i, Params[0], Params[1], Params[2], Params[3]))
                    i += 1

            if self.Num_Impropers > 0:
                File.write('\n\nImproper Coeffs # harmonic\n\n')
                i = 1
                for Params in self.Improper_Params:
                    File.write('%d %.4f -1 2\n' % (i, Params[0]))
                    i += 1

            File.write('\n\nAtoms # full\n\n')
            i = 1
            j = 1
            for Molecule_Obj in self.Molecule_List:
                for Atom_Obj in Molecule_Obj.Atom_List:
                    Atom_Obj.System_ID = i
                    File.write('%d %d %d %.8f %.4f %.4f %.4f\n' % (Atom_Obj.System_ID, Molecule_Obj.Mol_ID, Atom_Obj.LAMMPS_Type, Atom_Obj.Charge, Atom_Obj.Position[0], Atom_Obj.Position[1], Atom_Obj.Position[2]))
                    i += 1
                j += 1


            File.write('\n\nBonds\n\n')
            i = 1
            for Molecule_Obj in self.Molecule_List:
                for Bond_Obj in Molecule_Obj.Bond_List:
                    Bond_Obj.System_ID = i
                    
                    File.write('%d %d %d %d\n' % ( Bond_Obj.System_ID, Bond_Obj.LAMMPS_Type, Bond_Obj.Bond_Master.System_ID, Bond_Obj.Bond_Slave.System_ID))
                    i += 1


            File.write('\n\nAngles\n\n')
            i = 1
            for Molecule_Obj in self.Molecule_List:
                for Angle_Obj in Molecule_Obj.Angle_List:
                    Angle_Obj.System_ID = i
                    File.write('%d %d %d %d %d\n' % (Angle_Obj.System_ID, Angle_Obj.LAMMPS_Type, Angle_Obj.Angle_Slave1.System_ID, Angle_Obj.Angle_Master.System_ID,  Angle_Obj.Angle_Slave2.System_ID))
                    i += 1
            if self.Num_Dihedrals > 0:
                File.write('\n\nDihedrals\n\n')
                i = 1
                for Molecule_Obj in self.Molecule_List:
                    for Dihedral_Obj in Molecule_Obj.Dihedral_List:
                        Dihedral_Obj.System_ID = i
                        File.write('%d %d %d %d %d %d\n' % (Dihedral_Obj.System_ID, Dihedral_Obj.LAMMPS_Type, Dihedral_Obj.Dihedral_Slave1.System_ID, Dihedral_Obj.Dihedral_Master1.System_ID, Dihedral_Obj.Dihedral_Master2.System_ID, Dihedral_Obj.Dihedral_Slave2.System_ID))
                        i += 1

            if self.Num_Impropers > 0:
                File.write('\n\nImpropers\n\n')
                i = 1
                for Molecule_Obj in self.Molecule_List:
                    for Improper_Obj in Molecule_Obj.Improper_List:
                        Improper_Obj.System_ID = i
                        File.write('%d %d %d %d %d %d\n' % (Improper_Obj.System_ID, Improper_Obj.LAMMPS_Type, Improper_Obj.Improper_Master.System_ID, Improper_Obj.Improper_Slave1.System_ID, Improper_Obj.Improper_Slave2.System_ID, Improper_Obj.Improper_Slave3.System_ID))
                        i += 1
            # print "Moving lammp file to: " + self.ResultPath + self.Data_File
            # os.rename(self.Data_File, self.ResultPath + self.Data_File)
            # print self.ResultPath +"/" self.Data_File
        return

    # returns the dynamic type elements 1,2,3,...(gold_type-1)
    def dynamictype(self, gold_type):
        ret = ""
        for i in range(gold_type-1):
            ret += str(i+1)+" "
        return ret

    def Run_Lammps_Init(self, Nodes = 1):
        cmd = "mkdir " + Configure.Comet_Path % self.Name
        subprocess.call(["ssh", Configure.Comet_Login, cmd])
        
        # Set up input file
        In_Temp = Configure.Template_Path + "in.fixed_gold"
        print "In_Temp is: " + In_Temp
        In_File = "in.init_%s" % self.Name
        Sim_Name = "init_%s" % self.Name
        with open(In_Temp) as f:
            template = f.read()
        # print "Name is: " + self.Name

        # print "Atom params: " + str(self.Atom_Params)
        # print "Solvent log: " + str(self.Solvent_Log)
        # print "Mol List: " + str(self.Molecule_List)


        pair_string = ""
        i = 1
        for Params in self.Atom_Params:
            pair_string += 'pair_coeff %d %d lj/cut/coul/long %.3f %.3f\n' % (i, i, Params[2], Params[1])
            i += 1

        print "Pair string is: \n" + pair_string
        # print "atom_params: " + str(self.Atom_Params)
        gold_type = len(self.Atom_Params)
        print "gold type is: " + str(gold_type)
        dynamic_type = self.dynamictype(gold_type)
        s = template.format(System_Name = self.Name, pair_string = pair_string, gold_type = gold_type, dynamic_type = dynamic_type)
        with open(In_File,'w') as f:
            f.write(s)
        """
        Current Rule of thumb for Parallel Job Submission:
        MPI: 1 processor per 1000 atoms
        MPI + GPU: 4 GPU and 24 Processors --> Only use for >100K particles
        Need to do more benchmarks with current system --> good task for CJ
        """
        sub_temp = Configure.Template_Path + "sub_Lammps"
        NProcs = 12
        submit = "sub_%s" % self.Name
        with open(sub_temp) as f:
            template = f.read()
            s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(submit,'w') as f:
            f.write(s)
    
            
            
        File_Out1 = 'log.%s' % Sim_Name
        File_Out = 'restart.%s_init_1' % self.Name
        
        # Copy over to Comet
        os.system( Configure.c2c % (submit, self.Name))
        os.system( Configure.c2c % (In_File, self.Name))
        os.system( Configure.c2c % (self.Data_File, self.Name))
        os.system( Configure.c2l % (self.Name, File_Out1))
        
        try:
            File = open(File_Out1,'r')
        except:
            subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])

        Finished = False
        i = 0
        while not Finished:
            os.system( Configure.c2l % (self.Name, File_Out))
            try:
                File = open(File_Out,'r')
                Finished = True
            except:
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10
        os.system( 'rm %s' % File_Out)
        self.Current_Restart = File_Out
        return






    def Run_Lammps_NPT(self, GPU = False, Temp_Out = 0.0, time_steps = 1500000, Nodes = 1):
        """
            Function for running NPT dynamics for the system in lammps
        
        """
        Temp_In = self.Temperature
        count = int(self.Current_Restart.split('_')[-1])
        count += 1
        NPT_Temp = Configure.Template_Path + "in.SAM_NPT_Temp"
        NPT_In = "in.NPT_%s_%d_%d" % (self.Name, count, Temp_In)
        Sim_Name = NPT_In.split('.')[-1]
        if Temp_Out != 0.0:
            NPT_In += "_%d" % Temp_Out
            Sim_Name = NPT_In.split('.')[-1]
            self.Temperature = Temp_Out
            New_Restart = 'restart.%s_%d_%d' % (self.Name, Temp_Out, count)
        if Temp_Out == 0.0:
            Temp_Out = Temp_In
            New_Restart = 'restart.%s_%d_%d' % (self.Name, Temp_Out, count)
        with open(NPT_Temp) as f:
                template = f.read()
        pair_string = ""
        i = 1
        for Params in self.Atom_Params:
            pair_string += 'pair_coeff %d %d lj/cut/coul/long %.3f %.3f\n' % (i, i, Params[2], Params[1])
            i += 1
        print "atom_params: " + str(self.Atom_Params)
        gold_type = len(self.Atom_Params)
        print "gold type is: " + str(gold_type)
        dynamic_type = self.dynamictype(gold_type)
        s = template.format(Name = self.Name, pair_string = pair_string, gold_type = gold_type, dynamic_type=dynamic_type, count = count, Temp_In = Temp_In, Temp_Out = Temp_Out, Restart_In = self.Current_Restart, Restart_Out = New_Restart, Steps = time_steps)
        with open(NPT_In,'w') as f:
            f.write(s)
        sub_temp = Configure.Template_Path + "sub_Lammps"
        NProcs = 12
        submit = "sub_%s" % self.Name
        with open(sub_temp) as f:
            template = f.read()
            s = template.format(Sim_Name = Sim_Name, path = Configure.Comet_Path % self.Name, NProcs = Nodes*NProcs, Nodes=Nodes, tpn = NProcs)
        with open(submit,'w') as f:
            f.write(s)
            
        File_Out1 = 'log.%s' % Sim_Name
        File_Out = New_Restart
        Traj_File = self.Name + "_%d.lammpstrj" % count
            
        # Copy over to Comet
        os.system( Configure.c2c % (submit, self.Name))
        os.system( Configure.c2c % (NPT_In, self.Name))
        os.system( Configure.c2l % (self.Name, File_Out1))
            
        try:
            File = open(File_Out1,'r')
        except:
            subprocess.call(["ssh", Configure.Comet_Login, Configure.SBATCH % (self.Name, submit)])
            
        Finished = False
        i = 0
        while not Finished:
            os.system( Configure.c2l % (self.Name, File_Out))
            try:
                File = open(File_Out,'r')
                Finished = True
            except:
                print "Sleeping process", i, "minutes"
                time.sleep(600)
                i += 10
        
        os.system( 'rm %s' % File_Out)
        self.Current_Restart = File_Out
        os.system( Configure.c2l % (self.Name,  Traj_File))
        return


