#! usr/bin/python

# Import relevant modules
import numpy as np
import os
import subprocess
import shutil
import time
import random

# import Class structure
import Atom
import Bond
import Angle
import Dihedral
import Improper

Element_Dict = { 12.01078:"C", 1.0079469999999999:"H", 18.998:"F", 15.999:"O", 15.99943:"O", 32.06:"S", 14.007:"N" , 28.085529999999999:"Si", 196.97:"Au",196.96:"Au"}

class Surface(object):
    """
        Class defining a Surface produced from Andrew Summer's code
        instance variables (self) :
        N = Number of Atoms
        Name = Surface Name
        Atom_list[N] = List of Atom objects
        Bond_List[Num_Bonds] = list of bond objects
        Angle_List[Num_Angles] = list of angle objects
        Dihedral_List[Num_Dihedrals] = list of Dihedral objects
        Improper_List[Num_Impropers] = list of Improper objects
        self.Mol_ID = id for molecule reference
        self.Box_Bounds = bounds of simulation box
        """
    def __init__(self, Filename, Anneal=False):
        
        
        self.Name = Filename.split('.')[1].split('_')[0]
        File = open(Filename,'r')
        File_Lines = File.readlines()
        self.Bond_List = []
        self.Angle_List = []
        self.Dihedral_List = []
        self.Improper_List = []
        self.Mol_ID =0
        self.MW = 0.0
        self.basis = np.eye(3)
        self.Box_Length = []
        self.Molecule = False

        j = -1
        for Line in File_Lines:
            j += 1
            Line = Line.strip('\n').split()

            if len(Line) > 1:
                if Line[1] == 'atoms':
                    self.N = int(Line[0])
                    print self.N, "Atoms"
                if Line[1] == 'atom':
                    Atom_Types = int(Line[0])
                    self.Atom_List = np.empty(self.N, dtype=object)
                    Mass_List = np.empty(Atom_Types, dtype=float)
                    Pair_Coeffs = np.empty((Atom_Types,2), dtype=float)
                if Line[1] == 'bonds':
                    Num_Bonds = Line[0]
                if Line[1] == 'bond':
                    Bond_Types = int(Line[0])
                    Bond_Coeffs = np.empty((Bond_Types,2), dtype=float)
                if Line[1] == 'angles':
                    Num_Angles = Line[0]
                if Line[1] == 'angle':
                    Angle_Types = int(Line[0])
                    Angle_Coeffs = np.empty((Angle_Types,2),dtype=float)
                if Line[1] == 'dihedrals':
                    Num_Dihedrals = Line[0]
                if Line[1] == 'dihedral':
                    Dihedral_Types = int(Line[0])
                    Dihedral_Coeffs = np.empty((Dihedral_Types,4),dtype=float)
                if Line[1] == 'impropers':
                    Num_Impropers = Line[0]
                if Line[1] == 'improper':
                    Improper_Types = int(Line[0])
                    Improper_Coeffs = np.empty((Improper_Types,3),dtype=float)
                try:
                    if Line[2] == 'xlo':
                        self.Box_Length.append(float(Line[1]))
                    if Line[2] == 'ylo':
                        self.Box_Length.append(float(Line[1]))
                        self.Box_Length.append(100.0)
                except:
                    continue
                

                if Line[0] == 'Pair':
                    print "Getting Pair Coefficients"
                    for i in range(Atom_Types):
                        Pair_Coeffs[i,0] = float(File_Lines[i+j+2].split('\t')[1])
                        Pair_Coeffs[i,1] = float(File_Lines[i+j+2].split('\t')[2])
                if Line[0] == 'Bond':
                    print "Getting Bond Coefficients"
                    for i in range(Bond_Types):
                        Bond_Coeffs[i,0] = float(File_Lines[i+j+2].split('\t')[1])
                        Bond_Coeffs[i,1] = float(File_Lines[i+j+2].split('\t')[2])
                if Line[0] == 'Angle':
                    print "Getting Angle Coefficients"
                    for i in range(Angle_Types):
                        Angle_Coeffs[i,0] = float(File_Lines[i+j+2].split('\t')[1])
                        Angle_Coeffs[i,1] = float(File_Lines[i+j+2].split('\t')[2])
                if Line[0] == 'Dihedral':
                    print "Getting Dihedral Coefficients"
                    for i in range(Dihedral_Types):
                        Dihedral_Coeffs[i,0] = float(File_Lines[i+j+2].split('\t')[1])
                        Dihedral_Coeffs[i,1] = float(File_Lines[i+j+2].split('\t')[2])
                        Dihedral_Coeffs[i,2] = float(File_Lines[i+j+2].split('\t')[3])
                        Dihedral_Coeffs[i,3] = float(File_Lines[i+j+2].split('\t')[4])


            if len(Line) == 1:
                print Line
                if Line == ['']:
                    continue
                if Line == ['Atoms']:
                    print 'Getting Atoms'
                    for i in range(self.N):
                        Temp_Position = [float(File_Lines[i+j+2].split()[2]), float(File_Lines[i+j+2].split()[3]), float(File_Lines[i+j+2].split()[4]) ]
                        I_Flags = [0,0,0]
                        Temp_Position = np.asarray(Temp_Position, dtype=float)
                        
                        I_Flags = np.asarray(I_Flags, dtype=int)
                        self.Box_Length = np.asarray(self.Box_Length, dtype=float)
                        #print "IMAGE FLAGS", I_Flags
                        Temp_Position += I_Flags*self.Box_Length
                        Temp_Charge = 0.0
                        Type = int(float(File_Lines[i+j+2].split()[1]))
                        Mass = Mass_List[Type-1]
                        Element = Element_Dict[Mass]
                        Atom_Id = int(File_Lines[i+j+2].split()[0])
                        #print Element, Mass, Atom_Id
                        self.Atom_List[i] = Atom.Atom(Temp_Position, Element, Atom_Id)
                        self.Atom_List[i].Charge = Temp_Charge
                        self.Atom_List[i].Epsilon = Pair_Coeffs[Type-1,0]
                        self.Atom_List[i].Sigma = Pair_Coeffs[Type-1,1]
                        self.Atom_List[i].OPLS_Type = 1000 + Type
                    
                    self.Atom_List = sorted(self.Atom_List, key=lambda AtomO: AtomO.Atom_ID)
                    

                if Line == ['Masses']:
                    print "Extracting Masses"
                    for i in range(Atom_Types):
                        Mass_List[i] = float(File_Lines[i+j+2].split()[1])
                if Line == ['Bonds']:
                    print "Extracting Bonds"
                    for i in range(int(Num_Bonds)):
                        Bond_Info = File_Lines[i+j+2].strip('\n').split('\t')
                        Master = self.Atom_List[int(Bond_Info[2])-1]
                        Slave = self.Atom_List[int(Bond_Info[3])-1]
                        Kb = Bond_Coeffs[int(Bond_Info[1])-1, 0]
                        Req = Bond_Coeffs[int(Bond_Info[1])-1,1]
                        Bond_ID = int(Bond_Info[0])
                        self.Bond_List.append(Bond.Bond(Master, Slave,Req))
                        self.Bond_List[i].kb = Kb
                        self.Bond_List[i].Bond_ID = Bond_ID
                if Line == ['Angles']:
                    print "Extracting Angles"
                    for i in range(int(Num_Angles)):
                        Angle_Info = File_Lines[i+j+2].strip('\n').split('\t')
                        Slave1 = self.Atom_List[int(Angle_Info[2])-1]
                        Master = self.Atom_List[int(Angle_Info[3])-1]
                        Slave2 = self.Atom_List[int(Angle_Info[4])-1]
                        Ka = Angle_Coeffs[int(Angle_Info[1])-1,0]
                        Th0 = Angle_Coeffs[int(Angle_Info[1])-1,1]
                        Angle_ID = int(Angle_Info[0])
                        self.Angle_List.append(Angle.Angle(Master, Slave1, Slave2, Th0))
                        self.Angle_List[i].ka = Ka
                        self.Angle_List[i].Angle_ID = Angle_ID

                if Line == ['Dihedrals']:
                    print "Extracting Dihedrals"
                    for i in range(int(Num_Dihedrals)):
                        Dihedral_Info = File_Lines[i+j+2].strip('\n').split('\t')
                        Slave1 = self.Atom_List[int(Dihedral_Info[2])-1]
                        Master1 = self.Atom_List[int(Dihedral_Info[3])-1]
                        Master2 = self.Atom_List[int(Dihedral_Info[4])-1]
                        Slave2 = self.Atom_List[int(Dihedral_Info[5])-1]
                        Coeffs = Dihedral_Coeffs[int(Dihedral_Info[1])-1]
                        Dihedral_ID = int(Dihedral_Info[0])
                        self.Dihedral_List.append(Dihedral.Dihedral(Master1, Master2, Slave1, Slave2, 0.0))
                        self.Dihedral_List[i].Coeffs = Coeffs
                        self.Dihedral_List[i].Dihedral_ID = Dihedral_ID

                if Line == ['Impropers']:
                    print "Extracting Impropers"
                    for i in range(int(Num_Impropers)):
                        Improper_Info = File_Lines[i+j+2].strip('\n').split('\t')
                        Master = self.Atom_List[int(Improper_Info[2])-1]
                        Slave1 = self.Atom_List[int(Improper_Info[3])-1]
                        Slave2 = self.Atom_List[int(Improper_Info[4])-1]
                        Slave3 = self.Atom_List[int(Improper_Info[5])-1]
                        Coeff = Improper_Coeffs[int(Improper_Info[1])-1,0]
                        Improper_ID = int(Improper_Info[0])
                        self.Improper_List.append(Improper.Improper(Master, Slave1, Slave2, Slave3, Coeff, 180.0, Improper_ID))
                        
                        
        """
        # Compute center of Mass
        Mass_Weighted_Sum = np.zeros(3,dtype=float)
        for Atom_Obj in self.Atom_List:
            Mass_Weighted_Sum += Atom_Obj.Position*Atom_Obj.Mass
            self.MW += Atom_Obj.Mass
        
        print "Molecular Weight is ", self.MW, "Grams/Mole"
        self.COM = Mass_Weighted_Sum/self.MW
        print "COM is", self.COM
        
        # Zero COM
        for Atom_Obj in self.Atom_List:
            Atom_Obj.Position -= self.COM
        
        self.COM -= self.COM
        print self.COM
        """
        for Atom_Obj in self.Atom_List:
            self.MW += Atom_Obj.Mass
        return
        




