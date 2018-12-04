#! usr/bin/python

# OPLS

import Molecule
import Atom
import Configure

# Open up OPLS FILE
OPLS_Path = Configure.Template_Path + "oplsaa.prm.txt"
OPLS_FILE = open(OPLS_Path, 'r')
File_Lines = OPLS_FILE.readlines()

# Initialize Parameters Lists
OPLS_Bonds = []
OPLS_Angles = []
OPLS_Dihedrals = []
OPLS_Impropers = []
OPLS_VDW = []
OPLS_CHARGE = []

# Fill in parameter list
for Line in File_Lines:
    Line = Line.split()
    try:
        if Line[0] == "bond":
            OPLS_Bonds.append(Line)
        if Line[0] == "angle":
            OPLS_Angles.append(Line)
        if Line[0] == "torsion":
            OPLS_Dihedrals.append(Line)
        if Line[0] == "imptors":
            OPLS_Impropers.append(Line)
        if Line[0] == "vdw":
            OPLS_VDW.append(Line)
        if Line[0] == "charge":
            OPLS_CHARGE.append(Line)
    except:
        continue

def find_partial_charges(mol):
	for Atom_Obj in Molecule.Atom_List:
		for CHARGE in OPLS_CHARGE:
			if int(CHARGE[1]) == Atom_Obj.OPLS_Type:
				print Atom_Obj.OPLS_Type
				print "CHelpG:", Atom_Obj.Charge
				print "OPLS:", float(CHARGE[2])
				Atom_Obj.Charge = float(CHARGE[2])
