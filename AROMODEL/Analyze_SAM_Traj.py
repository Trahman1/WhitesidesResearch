#!usr/bin/python
# Import relevant modules
import numpy as np
# from numpy import lin
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
from DistanceMatrix import DistanceMatrix
import OPLS
import Atom as AtomOld
import Molecule as MoleculeOld
import DataFilesHandler

# TODO: Automate Mass Type
Masses = [32.060, 12.011, 12.011, 15.99, 12.011, 1.008, 196.960]
Elements = ['S', 'C', 'C', 'O', 'C', 'H', 'Au']

# Masses = [32.060, 12.011, 12.011, 12.011, 15.99, 14,007, 12.011, 12.011, 1.008, 196.960]

# The Amide_11_1 gives a problem here so if the masses are actually needed for something, that will need to be taken into account
rcParams['figure.figsize'] = 6, 5
mpl.rcParams['axes.labelsize'] = 'large'
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
        self.Element = Elements[Type-1]
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
        self.weightedPosition = np.empty(3, dtype = float)
        self.Atom_List = []
        self.dipoleMoment = np.empty(3, dtype = float)
        # self.COM = np.empty(3, dtype = float)
        return
    
    def Add_Atom(self, Atom):
        self.N += 1
        self.Atom_List.append(Atom)
        self.MW += Masses[Atom.Type-1]
        self.weightedPosition = Masses[Atom.Type-1] * Atom.position
        return
    
    def NP_Convert(self):
        self.Atom_List = np.asarray(self.Atom_List)
        return

    def COM(self):
        return self.weightedPosition/self.MW



class Snap_Shot(object):
    """
        Class defining a snap shot of the trajectory
        instance variables: Time_step, Mol_List, Box_Dim, Rg_Dist, E2E_Dist, RDF_Dist
    """
    def __init__(self, Time_Step, Mol_List, Box_Dim, Atom_List, N):
        self.Time_Step = Time_Step
        self.Mol_List = Mol_List
        self.Box_Dim = Box_Dim
        self.num_Mols = len(Mol_List)
        self.Atom_List = Atom_List
        self.DistanceMatrix = DistanceMatrix()
        self.DistanceMatrix.setPeriodicBound(Box_Dim)
        self.numCarbon = 5680.
        boxx, boxy, boxz = Box_Dim
        return
    
    def showDistanceDistribution(self):
        self.DistanceMatrix.showDistanceDistribution()

    def pairCorrelation_Hist(self, numbins=200):
        # print "Calculating Snapshot Distance Matrix"
        for Mol in self.Mol_List:
            self.DistanceMatrix.addPoint(Mol.COM())
        # print "Finished placing atoms...Calculating Distance Distribution"
        self.DistanceMatrix.calculateDistanceDistribution()
        # print "returning pyplot"
        n, bins = self.DistanceMatrix.pyplotDistanceDistribution(numbins = 200)
        print "n is: "
        print n
        print bins
        return n, bins

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

    def Compute_Tangent_Correlation(self):
        TangentCorrelation = []
        for Mol in self.Mol_List:
            MolTangentCorrelation = []
            RelevantAtoms = [Atom for Atom in Mol.Atom_List.tolist() if ((Atom.Mass > 2.) and (Atom.Mass < 196))]
            # print "mass of second atom is: " + str(Mol.Atom_List[1].Mass)
            initialvector = RelevantAtoms[1].position-RelevantAtoms[0].position
            initialvector /= np.linalg.norm(initialvector)
            lastpos = RelevantAtoms[0].position
            for Atom in RelevantAtoms[1:]:
                # print "adding new pos: " + str(Atom.position)
                newpos = Atom.position
                newvector = newpos-lastpos
                newvector /= np.linalg.norm(newvector)
                # print "tangent correlation is: " + str(self.Tangent_Correlation(latestpos[0],latestpos[1],newpos))
                MolTangentCorrelation.append(np.dot(initialvector,newvector))
                lastpos =newpos
            # print "Mol Tangent Correlation: " + str(np.asarray(MolTangentCorrelation)) + " with mass"
            TangentCorrelation.append(np.asarray(MolTangentCorrelation))
        TangentCorrelation = np.asarray(TangentCorrelation) 
        # print "TangentCorrelation is: " + str(TangentCorrelation)
        TangentCorrelation = np.mean(TangentCorrelation, axis = 0)
        # print "TangentCorrelation is: " + str(TangentCorrelation)
        return TangentCorrelation

    def carbonDensity_Hist(self, numbins = 200):
        Z_pos = []
        Z_max = []
        theta = []
        N = 0
        for Mol in self.Mol_List:
            mol_z = []
            for Atom in Mol.Atom_List:
                if Atom.Mass == 12.011:
                    # print "found Carbon"
                    N += 1
                    Z_pos.append(Atom.position[2])
                    mol_z.append(Atom.position[2])
            try:
                Z_max.append(np.asarray(mol_z).max())
            except:
                continue
        # Z_max_av = np.asarray(Z_max).mean()
        # Z_max_std = np.asarray(Z_max).std()
        n, bins = np.histogram(Z_pos, bins= numbins, range=(0.07,0.3))
        return n/self.numCarbon, bins*100-7.6
        # return Zn, Zbins, Z_max_av, Z_max_std

    def get_net_dipole(self, molclass, molName):
        # partial_charges = OPLS.get_partial_charges(molclass, molName)
        partial_charges = [(0,.32), (0,-.64)]
        print partial_charges
        net_dipole = np.zeros(3)
        # print pc
        # splicing out gold surface
        for mol in self.Mol_List:
            for atom in mol.Atom_List:
                if atom.Element == "Au":
                    continue
                else:
                    partial_charge = partial_charges[atom.Type-1][1]
                    # print atom.position
                    # print partial_charge
                    dipole_moment = partial_charge*atom.position
                    # print dipole_moment
                    net_dipole = net_dipole+dipole_moment
                    mol.dipoleMoment = mol.dipoleMoment+dipole_moment
        print "Net dipole is: ..." + str(net_dipole)
        # return net_dipole
        # for mol in self.Mol_List[:-1]:
        #     if (np.dot(mol.dipoleMoment,net_dipole)/(np.linalg.norm(mol.dipoleMoment)*np.linalg.norm(net_dipole)) < -.9):
        #         print "alert!"
        #         print mol.dipoleMoment
        #         for atom in mol.Atom_List:
        #             print atom.Element, str(atom.position)
        #         # print net_dipole
        #         break
        plt.hist([np.dot(mol.dipoleMoment,net_dipole)/(np.linalg.norm(mol.dipoleMoment)*np.linalg.norm(net_dipole)) for mol in self.Mol_List[:-1]], bins=100)
        print os.getcwd()
        plt.savefig(molclass+"_"+molName+"_dipoleDistribution.png")
        print "saved!"
        plt.clf()
        return net_dipole

                    # print atom.Type, atom.Element
                # else:
                #     print atom.Type, atom.Element
            # for atom in mol.Atom_List:
            #     if atom.Element == "Au":
            #         print "found: ", atom.Mol_ID
            #     else:
            #         print "not found", atom.Mol_ID
            #         # break




        # for Mol in self.Mol_List[:1]:
        #     # for Atom in Mol.Atom_List[:1]:
        #     #     print Atom.Type
        #     #     print Atom.Mass
        #     #     print Atom.id
        #     alist = []
        #     for index, atom in enumerate(Mol.Atom_List):
        #         alist.append(AtomOld.Atom(atom.position, atom.Element, index))
        #     m = MoleculeOld.Molecule("Test", None, alist)
        #     m.Set_Up_FF()
        #     OPLS.Assign_OPLS(m)

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
    def __init__(self, File_Name, molclass, molName, label=None):
        # label defaults to molName if not otherwise specified
        self.label = label
        if label == None:
            self.label = molName
        self.mol_class = molclass
        self.mol = molName

        self.SnapShotList = self.readSnapshots(File_Name)

        self.Num_Snaps = len(self.SnapShotList)
        self.graphsPath = Configure.molGraphsPath(molclass, molName)
        if not os.path.exists(self.graphsPath):
            os.mkdir(self.graphsPath)

        # TODO: Check numCarbons
        self.numCarbon = 5670
        return
    def changelabel(self, newlabel):
        self.label = newlabel

    def saveCarbonDensity(self):
        self.averageHistograms(lambda ss, numbins: ss.carbonDensity_Hist(numbins), "Carbon Density", numbins = 200,xlabel = "Number of Carbons", ylabel = "Thickness ($\AA$)")
        self.save_and_clear("CarbonDensity")

    def savePairCorrelation(self, label=None):
        ssmethod = lambda ss, numbinsprime : ss.pairCorrelation_Hist(numbins=numbinsprime)
        self.averageHistograms(ssmethod, "Pair Correlation", xlabel = "r", ylabel = "g(r)")
        plt.xlim(0,min(self.SnapShotList[0].Box_Dim[:2])/2)
        plt.ylim(bottom = .75)
        plt.hlines(1,0,min(self.SnapShotList[0].Box_Dim[:2])/2, colors = "gray", linestyles = "dashed")
        self.save_and_clear("PairCorrelation")

    def saveTangentCorrelation(self):
        self.averagePlots(lambda ss: ss.Compute_Tangent_Correlation(), "Tangent Correlation", xlabel = "Atom #")
        self.save_and_clear("TangentCorrelation")
        # print "TANGENT CORRELATION IS: " + str(totalTangentCorrelation)
        # plt.plot(totalTangentCorrelation, label = self.label)

    def save_and_clear(self, title):
        plt.savefig(Configure.molGraphsPath(self.mol_class, self.mol) + self.mol + "_" + title)
        plt.clf()

    def averageHistograms(self, snapShotMethod, title, numbins = 200, label=None, xlabel = None, ylabel= None):
        if label is None:
            label = self.label
        n, __bins, __patches = 0, None, None
        for ss in self.SnapShotList:
            ntemp, bintemp = snapShotMethod(ss, numbins)
            n = n+ntemp
        plt.plot(bintemp[0:-1], [x/self.Num_Snaps for x in n], label = label)
        plt.legend(loc = 'upper right', frameon = False, fontsize= 30)
        plt.title(label + " " + title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

    def averagePlots(self, snapShotMethod, title, label=None, xlabel = None, ylabel = None):
        total = []
        for ss in self.SnapShotList:
            total.append(np.asarray(snapShotMethod(ss)))
        total = np.asarray(total)
        total = total.mean(axis=0)
        if (label == None):
            label = self.label
        plt.plot(total, label = label)
        plt.legend( loc = 'upper right', frameon = False, fontsize= 30)
        plt.title(label + " " + title)
        plt.xlabel(xlabel)

    

    # parses through file and creates snapshot list
    def readSnapshots(self, File_Name):

        with open(File_Name) as f:
            File_Lines = f.readlines()
        Temp_Snap_List = []
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
                    # print 
                    Atom_Line = File_Lines[i + 9 + j].split(' ')
                    ID = int(Atom_Line[0])
                    TYPE = int(Atom_Line[1])
                    MOL = int(Atom_Line[2])
                    POS= [float(Atom_Line[3]), float(Atom_Line[4]), float(Atom_Line[5])]
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
                Snap_Shot_Obj = Snap_Shot(Time_Step, Mol_List, Box_Dim, Sorted_Atom_List, N)
                Temp_Snap_List.append(Snap_Shot_Obj)
        
        return Temp_Snap_List  

def savePlots(molclass, mol, traj):
    newtraj = get_trajectory_file(molclass, mol, traj)
    newtraj.saveCarbonDensity()
    newtraj.saveTangentCorrelation()
    newtraj.savePairCorrelation()

def get_trajectory_file(molclass, mol, traj):
    print "reading trajectory for: " + molclass + ": " + mol + ": " + str(traj)
    respath = Configure.Results_Path + molclass + "/"
    molrespath = respath + mol + "/"
    trajpath = molrespath+mol+"_" + str(traj) + ".lammpstrj"
    if (traj == "1"):
        print "traj adjusted"
        trajpath = molrespath+"data."+mol
    newtraj = Trajectory(trajpath,molclass,mol, label=mol)
    return newtraj

def main():
    scriptname, molclass, mol, trajnum = sys.argv 
    
    # test()

    # files = [("OEG1", 6),("OEG2", 6),("OEG3", 6),("OEG4", 6),("OEG5", 6),("OEG6", 6)]
    # for mol, trajnum in files:
    savePlots(molclass, mol, trajnum)
    # traj.savePlots()
    #     ss = traj.SnapShotList[0]
    #     ss.get_net_dipole(molclass, mol)


def test():
    H2Otest()

def H2Otest():
    H1 = Atom(0, 1, 1, np.asarray([1.81,0,0]), np.asarray([0,0,0]))
    H2 = Atom(1, 1, 1, np.asarray([0,1.81,0]), np.asarray([0,0,0]))
    O = Atom(2, 2, 1, np.asarray([0,0,0]), np.asarray([0,0,0]))
    h2O = Molecule(1)
    h2O.Add_Atom(H1)
    h2O.Add_Atom(H2)
    h2O.Add_Atom(O)
    Atom_List = h2O.Atom_List
    Mol_List = np.asarray([h2O])
    boxDim = np.asarray([10,10,10])
    ss = Snap_Shot(0, Mol_List, boxDim, Atom_List, 3)
    ss.get_net_dipole("Test", "H2O")

    # savePlots(molclass, mol, trajnum)

    # for traj in trajnums:
    #     newtraj = Trajectory(molrespath+mol+"_" + str(traj) + ".lammpstrj",molclass,mol, label=mol)
    #     newtraj.savePlots()
        # newtraj.Compute_Tangent_Correlation()
    # plt.show()
    # for Traj in TrajObjects:

    # plt.legend( loc = 'upper right', frameon = False, fontsize= 20)
    # plt.show()
    # plotname = "TangentCorrelation"
    # plt.savefig(Configure.Results_Path + "/" + molclass + "/" + molclass + plotname)

    # plt.savefig("molclassCarbonDensity.png")
    
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



