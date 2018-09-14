#!usr/bin/python
# Import relevant modules
import numpy as np
import matplotlib.pyplot as plt
import math
import pickle
import matplotlib as mpl
from pylab import rcParams
import time
from scipy.signal import argrelextrema
import sys
import Configure
import os

Masses = [32.060, 12.011, 12.011, 15.99, 12.011, 1.008, 196.960]
# The Amide_11_1 gives a problem here so if the masses are actually needed for something, that will need to be taken into account
rcParams['figure.figsize'] = 12, 10
mpl.rcParams['axes.color_cycle'] = ['b', 'k','r', 'c', 'g', 'm']

# Define Classes
class Atom(object):
    """
        class defining an atom
        instance variables: Type, id,  Mol_ID, position[3], image_flags[3]
        """
    def __init__(self, id, Type, Mol_ID, position, image_flags):
        self.Type = Type
        self.id = id
        self.position = np.asarray(position, dtype= float)
        self.image_flags = np.asarray(image_flags, dtype= int)
        self.Mol_ID = Mol_ID
        self.unwrapped_position = np.zeros(3, dtype =float)
        self.Mass = Masses[Type-1]
        return


    def Print_Info(self):
        print "Atom ID = %d, Type = %d\n, Mol_ID = %d" % (self.id, self.Type, self.Mol_ID)
        print "Position"
        print self.position
        print self.image_flags
        print "-----------------------------------------"
        return


class Molecule(object):
    """
        class defining a molecule
        instance variable: Mol_ID, N, MW, COM, RG
        Note: COM and RG are computed as mass averages
        """
    def __init__(self, Mol_ID):
        self.Mol_ID = Mol_ID
        self.N = 0
        self.MW = 0.0
        self.COM = np.empty(3, dtype = float)
        self.Atom_List = []
        return
    
    def Add_Atom(self, Atom):
        self.N += 1
        self.Atom_List.append(Atom)
        self.MW += Masses[Atom.Type-1]
        return
    
    def NP_Convert(self):
        self.Atom_List = np.asarray(self.Atom_List)
        return

class Snap_Shot(object):
    """
        Class defining a snap shot of the trajectory
        instance variables: Time_step, Mol_List, Box_Dim, Rg_Dist, E2E_Dist, RDF_Dist
    """
    def __init__(self, Time_Step, Mol_List, Box_Dim):
        self.Time_Step = Time_Step
        self.Mol_List = Mol_List
        self.Box_Dim = Box_Dim
        return
    
    def Compute_COM(self):
        """
            Compute the centers of mass of all the polymer chains
            """
        for Mol in self.Mol_List:
            Mass_Weighted_Sum = np.zeros(3,dtype=float)
            for Atom in Mol.Atom_List:
                Atom.unwrapped_position = Atom.position + np.multiply(self.Box_Dim, Atom.image_flags)
                Mass_Weighted_Sum += Polymer.Mass[Atom.Type]*Atom.unwrapped_position
            Mol.COM = Mass_Weighted_Sum/Mol.MW
        return

    def Compute_Histogram(self):
        Z_pos = []
        Z_max = []
        theta = []
        N = 0
        for Mol in self.Mol_List:
            mol_z = []
            for Atom in Mol.Atom_List:
                if Atom.Mass == 12.011:
                    N += 1
                    Z_pos.append(Atom.position[2])
                    mol_z.append(Atom.position[2])
            try:
                Z_max.append(np.asarray(mol_z).max())
            except:
                continue
                
        Z_max_av = np.asarray(Z_max).mean()
        Z_max_std = np.asarray(Z_max).std()
        n, bins = np.histogram(Z_pos, bins = 100, range=(0.07,0.3))
        #print n/float(N), bins
        #plt.hist(Z_pos, bins ='auto', normed=True,histtype='bar')
        #plt.legend( loc = 'upper right', frameon = False, fontsize= 25)
        #plt.tick_params( labelsize = 20, width=2, length=7)
        # plt.xlabel("Z height ($\AA$)", fontsize=25)
        #plt.ylabel("Probability", fontsize=25)
        #plt.show()
        return n, bins, Z_max_av, Z_max_std
"""
    def Compute_Angle(self):
        height = {}
        for Mol in self.Mol_List:
            for Atom in Mol.Atom_List:
"""


class Trajectory(object):
    """
        Class defining a complete trajectory outputted from an MD simulation
        instance variables: Num_Snaps, Snap_List
        """
    
    def __init__(self, File_Name, molName):
        File = open(File_Name,'r')
        File_Lines = File.readlines()
        Temp_Snap_List = []
        n = np.zeros(100, dtype=float)

        # Parse through file lines and extract the trajectory data
        for i in range(len(File_Lines)):
            line = File_Lines[i]
            if line == "ITEM: TIMESTEP\n":
                Time_Step = int(File_Lines[i+1])
                N = int(File_Lines[i+3])
                Atom_List = np.empty(N, dtype=object)
                XBounds = [ float(File_Lines[i+5].split()[0]), float(File_Lines[i+5].split()[1])]
                YBounds = [ float(File_Lines[i+6].split()[0]), float(File_Lines[i+6].split()[1])]
                ZBounds = [ float(File_Lines[i+7].split()[0]), float(File_Lines[i+7].split()[1])]
                Box_Dim = [ abs(XBounds[1]- XBounds[0]), abs(YBounds[1] - YBounds[0]), abs(ZBounds[1]- ZBounds[0])]
                
                #print "Timestep = %d, N = %d\n" % (Time_Step, N)
                #print Box_Dim
                # Extract Snapshot as list of atoms
                #print "Extracting Atoms...\n"
                for j in range(N):
                    Atom_Line = File_Lines[i + 9 + j].split(' ')
                    ID = int(Atom_Line[0])
                    TYPE = int(Atom_Line[1])
                    MOL = int(Atom_Line[2])
                    POS= [ float(Atom_Line[3]), float(Atom_Line[4]), float(Atom_Line[5])]
                    try:
                        IMAGE = [ int(Atom_Line[6]), int(Atom_Line[7]), int(Atom_Line[8])]
                    except:
                        IMAGE = [0,0,0]
                    Atom_List[j] = Atom(ID, TYPE, MOL, POS, IMAGE)
                #print "Finished extracting atoms, now sorting them"
                Sorted_Atom_List = sorted(Atom_List, key=lambda Atom: Atom.id)
                N_MOL = Sorted_Atom_List[-1].Mol_ID
                # Extract Molecule Objects from Atom List
                #print "Instantiating Molecule Objects"
                #print "N_MOL = %d " % N_MOL
                Mol_List = np.empty(N_MOL, dtype=object)
                for i in range(N_MOL):
                    Mol_List[i] = Molecule(i+1)
                
                
                for Atom_Obj in Sorted_Atom_List:
                    #Atom_Obj.Print_Info()
                    MOLID = Atom_Obj.Mol_ID
                    Mol_List[MOLID-1].Add_Atom(Atom_Obj)
                    Atom_Obj.position -= [ XBounds[0], YBounds[0], ZBounds[0]]
                #print "Converting to Numpy array"
                for Mol_Obj in Mol_List:
                    Mol_Obj.NP_Convert()
                
                #print "Instantiating Snap Shot Object, Computing properties"
                Snap_Shot_Obj = Snap_Shot(Time_Step, Mol_List, Box_Dim)
                ntemp, bins, z_max_av, z_max_std = Snap_Shot_Obj.Compute_Histogram()
                n += ntemp
                Temp_Snap_List.append(Snap_Shot_Obj)
        self.z_max_av = z_max_av*100 - 7.6
        self.z_max_std = z_max_std*100
        # Set instance variables
        self.Num_Snaps = len(Temp_Snap_List)
        self.Snap_Shot_List = Temp_Snap_List
        plt.plot(bins[0:-1]*100 - 7.6, n/(float(self.Num_Snaps)*5670.), linewidth=5, label=molName)
        for i in range(1,len(n)):
            if n[i] == 0.0 and n[i-1] != 0.0:
                self.height = (bins[i])*100 - 7.6
        print self.height
        extrema = argrelextrema(n, np.greater)
        self.height = bins[extrema[0][-1]]*100 - 7.6
        print extrema
        print self.height
        #plt.axvline(x=self.height, linestyle = '--', color ='k', alpha='0.5', linewidth=5)
        plt.legend( loc = 'upper right', frameon = False, fontsize= 30)
        plt.tick_params( labelsize = 30, width=2, length=7)
        plt.xlabel("Z height ($\AA$)", fontsize=30)
        plt.ylabel("Carbon Atom Density (1/$\AA^2$)", fontsize=30)
        plt.xlim((0,25))
        return



def main():
    scriptname, molclass = sys.argv 
    respath = Configure.Results_Path + molclass + "/"
    TrajObjects = []
    for i in range(1,8):
        mol = molclass + str(i)
        print "Adding mol: " + str(mol)
        molrespath = respath + mol + "/"
        latest_traj = 0
        # for j in range(1,7):
        #     if (os.path.exists(molrespath+mol+"_" + str(j) + ".lammpstrj")):
        #         latest_traj = j
        # print "Latest traj is: " + str(latest_traj)
        try:
            TrajObjects.append(Trajectory(molrespath+mol+"_3" + ".lammpstrj", mol))
            print "Success"
        except:
            print "Did not find file"
            continue
        
    plt.show()
    plt.savefig("molclassCarbonDensity.png")
    
    # numc = []
    # height = []
    # height_std = []
    # numlist = [[12,1], [10,1], [8,1], [6,1], [4,1]]
    # for i in numlist:
    #     #try:
    #     print 'Amide_%d_%d/Amide_%d_%d_5.lammpstrj' %(i[0],i[1],i[0],i[1])
    #     Traj_Object1 = Trajectory('Amide_%d_%d/Amide_%d_%d_5.lammpstrj' %(i[0],i[1],i[0],i[1]))
    #     height.append(Traj_Object1.z_max_av)
    #     height_std.append(Traj_Object1.z_max_std)
    #     numc.append(i[0])
    #     #except:
    #     #    print "couldn't find", i
    #     #   continue

    # plt.show()
    # height = np.asarray(height)
    # numc = np.asarray(numc)
    # numc2 = [6,8,10,12]
    # exp = [ 11, 13.9, 14.9, 17.1]
    # Yerr = [1.5,1.5,1.5,1.5]
    # plt.errorbar(numc, height, yerr=height_std, marker='o', markersize=20, linewidth=5, color='k', capsize=5,elinewidth=5,markeredgecolor = 'k', markeredgewidth=3, label='Simulated')
    # plt.errorbar(numc2, exp, yerr=Yerr, marker='o',linestyle='None', markersize=20, capsize=5,elinewidth=5,markeredgecolor = 'r', markeredgewidth=3, color='r', label='Measured' )

    # plt.legend( loc = 'upper left', frameon = False, fontsize= 30)
    # plt.xlim((2, 14))
    # plt.tick_params( labelsize = 30, width=2, length=7)
    # plt.xlabel("Number of Carbons", fontsize=30)
    # plt.ylabel("Thickness ($\AA$) ", fontsize=30)
    # print height
    # plt.show()

if __name__=='__main__': main()



