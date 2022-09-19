import numpy as np
from collections import Counter
from ase.data import covalent_radii as CR
from ce_expansion.atomgraph.bcm import BCModel
import ase.cluster as ac
from ase.visualize import view
import sys
import os
import argparse


# HELPER FUNCTIONS
def Dir_Exist_Or_Create(MYDIR):
    """If a directory exists then do nothing.  If not then create it.  
    Args:
        MYDIR (str): The name of the directory you are checking for
    """
    MYDIR = (MYDIR)
    CHECK_FOLDER = os.path.isdir(MYDIR)
    # If folder doesn't exist, then create it.
    if not CHECK_FOLDER:
        os.makedirs(MYDIR)
        current_directory = os.getcwd()
        print("created folder : ", os.path.join(current_directory,MYDIR))
        print('\n')

def get_coordination_numbers(atoms, covalent_percent=1.25):
    """Returns an array of coordination numbers and an array of existing bonds determined by
    distance and covalent radii.  By default a bond is defined as 120% of the combined radii
    or less. This can be changed by setting 'covalent_percent' to a float representing a 
    factor to multiple by (default = 1.2).

    If 'exclude' is set to an array,  these atomic numbers with be unable to form bonds.
    This only excludes them from being counted from other atoms,  the coordination
    numbers for these atoms will still be calculated,  but will be unable to form
    bonds to other excluded atomic numbers.
    """

    # Get all the distances
    distances = np.divide(atoms.get_all_distances(mic=True), covalent_percent)
    
    # Atomic Numbers
    numbers = atoms.numbers
    # Coordination Numbers for each atom
    cn = []
    cr = np.take(CR, numbers)
    # Array of indices of bonded atoms.  len(bonded[x]) == cn[x]
    bonded = []
    indices = list(range(len(atoms)))
    for i in indices:
        bondedi = []
        for ii in indices:
            # Skip if measuring the same atom
            if i == ii:
                continue
            if (cr[i] + cr[ii]) >= distances[i,ii]:
                bondedi.append(ii)
        # Add this atoms bonds to the bonded list
        bonded.append(bondedi)
    for i in bonded:
        cn.append(len(i))
    return cn, bonded

def evenly_distribute(atoms,CNs,CN,Atom_Type_2):
    """ Evenly distributes a second atom type over a specified coordination environment.

    Args:
        atoms (ase.Atoms): The atoms object to modify.
        CNs (dict): The output from get_coordination_numbers.
        CN (int): The coordination number to evenly distribute across.
        Atom_Type_2 (str): The atom type to add.

    Returns:
        atoms (ase.Atoms): The modified atoms object.
    """

    indices = [i for i, x in enumerate(CNs) if x == CN]
    for i,x in enumerate(indices):
        if i % 2:
            atoms[x].symbol = Atom_Type_2
    return atoms


if __name__ == "__main__":
    # Parse the command line arguments
    parser = argparse.ArgumentParser(description='Generates a set of NPs for gamma values')
    parser.add_argument('Atom_Type_1', type=str, help='The first atom type')
    parser.add_argument('Atom_Type_2', type=str, help='The second atom type')
    # make number of shells default to 4
    parser.add_argument('-n', '--number_of_shells', type=int, default=4, help='The number of shells to use default=4')

    args = parser.parse_args()
    Atom_Type_1 = args.Atom_Type_1
    Atom_Type_2 = args.Atom_Type_2
    number_of_shells = args.number_of_shells

    # CALCULATIONS
    Dir_Exist_Or_Create(f'Data')
    os.chdir('Data')
    folder_name = Atom_Type_2+Atom_Type_1
    Dir_Exist_Or_Create(folder_name)

    atoms = ac.Icosahedron(Atom_Type_1, number_of_shells)

    # Writing it to a file
    # atoms.write(os.path.join(folder_name,Atom_Type_1+str(len(atoms))+'.xyz'))


    CNs, Bonds = get_coordination_numbers(atoms)
    Unique_CNs = np.unique(CNs) # For icosahedron morphology

    for CN in Unique_CNs:
        atoms = evenly_distribute(atoms,CNs,CN,Atom_Type_2)

    ## Write to a file

    file_name_1 = Atom_Type_1+str(len(atoms[atoms.symbols==Atom_Type_1]))+Atom_Type_2+str(len(atoms[atoms.symbols==Atom_Type_2]))+'.xyz'
    atoms.set_initial_charges(CNs)
    atoms.write(os.path.join(folder_name,file_name_1))

    # Now we just swap the orderings

    indices_Atom_Type_1 = [i for i,x in enumerate(atoms.symbols) if x==Atom_Type_1]

    indices_Atom_Type_2 = [i for i,x in enumerate(atoms.symbols) if x==Atom_Type_2]

    atoms2 = atoms.copy()

    for i,x in enumerate(indices_Atom_Type_1):
        atoms2[x].symbol = Atom_Type_2
            
    for i,x in enumerate(indices_Atom_Type_2):
        atoms2[x].symbol = Atom_Type_1
    # Saving the second file
    atoms2.set_initial_charges(CNs)
    file_name_2 = Atom_Type_1+str(len(atoms2[atoms2.symbols==Atom_Type_1]))+Atom_Type_2+str(len(atoms2[atoms2.symbols==Atom_Type_2]))+'.xyz'
    atoms2.write(os.path.join(folder_name,file_name_2))
