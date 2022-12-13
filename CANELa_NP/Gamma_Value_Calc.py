import pandas as pd
import math
from ase.io import read
import numpy as np
from ce_expansion.atomgraph.bcm import BCModel
from BCM_Sandbox import BCM_Mod
from ce_expansion.atomgraph.adjacency import build_bonds_arr
from sympy.solvers import solve
from sympy import Symbol
import collections.abc
import json
from os.path import exists
import os
import argparse

# Constants 
"""Please add your own Bulk CE values here for whatever metals you are using"""
ce_bulk_pbe_d3 = {'Ag': 0.0, 'Al': 0.0, 'Au': 0.0, 'Cu': 0.0, 'Fe': 0.0, 'Ni': 0.0, 'Pd': 0.0, 'Pt': 0.0, 'Rh': 0.0, 'Ru': 0.0, 'Ti': 0.0, 'Zn': 0.0}
ce_bulk_pbe_d3 = {'Au':-3.64,'Pd':-4.20,'Pt':-6.20,"Ag":-2.96,"Cu":-3.95,"Ni":-5.11,"Ir":-7.95}
Ha_to_eV = 27.2114

# HELPER FUNCTIONS
def recursive_update(d: dict, u: dict) -> dict:
    """
    recursively updates 'dict of dicts'
    Ex)
    d = {0: {1: 2}}
    u = {0: {3: 4}, 8: 9}

    recursive_update(d, u) == {0: {1: 2, 3: 4}, 8: 9}

    Args:
    d (dict): the nested dict object to update
    u (dict): the nested dict that contains new key-value pairs

    Returns:
    d (dict): the final updated dict
    """
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = recursive_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def custom_calc_ce(atoms,gammas,ce_bulk):
    bcm_custom = BCM_Mod(atoms,gammas=gammas,ce_bulk=ce_bulk)
    return bcm_custom.calc_ce()

def calc_gammas(atoms,CEs):
    gamma_a = Symbol('gamma_a')
    gamma_b = Symbol('gamma_b')
    A,B = np.unique(atoms[0].symbols)
    gammas = {A:{A:1,B:gamma_a},B:{B:1,A:gamma_b}}

    sol = solve([custom_calc_ce(atoms[0],gammas,ce_bulk_pbe_d3) + custom_calc_ce(atoms[1],gammas,ce_bulk_pbe_d3) - sum(CEs),
                gamma_a+gamma_b - 2])
    # compare the custom calc ce to the ce values
    new_gammas = {A:{A:1,B:float(sol[gamma_a])},B:{B:1,A:float(sol[gamma_b])}}
    print("CEs: ",CEs)
    print("CEs from BCM: ",custom_calc_ce(atoms[0],new_gammas,ce_bulk_pbe_d3),custom_calc_ce(atoms[1],new_gammas,ce_bulk_pbe_d3))

    return new_gammas

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Calculate Gamma Values for a given set of CEs')
    # Parse the command line arguments
    parser.add_argument('Atom_Type_1', type=str, help='The first atom type')
    parser.add_argument('Atom_Type_2', type=str, help='The second atom type')
    parser.add_argument('-n', '--total_atoms', type=int, default=147, help='The total number of atoms = 147')
    
    # Collect the arguments
    args = parser.parse_args()
    Atom_Type_1 = args.Atom_Type_1
    Atom_Type_2 = args.Atom_Type_2
    number_of_atoms = args.total_atoms

    Atoms_folder = Atom_Type_2 + Atom_Type_1
    
    atoms_1_comp = Atom_Type_1 + str(math.floor(number_of_atoms/2)) + Atom_Type_2 + str(math.ceil(number_of_atoms/2))
    atoms_1_name = atoms_1_comp  + '.xyz'
    atoms_2_comp = Atom_Type_1 + str(math.ceil(number_of_atoms/2)) + Atom_Type_2 + str(math.floor(number_of_atoms/2)) # + '.xyz'
    atoms_2_name = atoms_2_comp + '.xyz'

    # Find paths to all the files
    gamma_folder_name = os.path.join(os.path.dirname(__file__), "Data") # Path to the folder containing all the data
    gamma_values_path = os.path.join(gamma_folder_name, "np_gammas.json") # Path to the file where the gamma values will be stored
    atom_energies = os.path.join(gamma_folder_name, "single_atom_energies.xlsx") # This should contain the single point energies of all the atoms in the system
    NP_energies = os.path.join(gamma_folder_name, "NP_energies.csv") # This should contain the energies of all the NP's you want to calculate the gamma values for
    NP1 = os.path.join(gamma_folder_name, Atoms_folder,atoms_1_name) # This should be the path to the A_{(n-1)/2}B_{(n+1)/2} NP file
    NP2 = os.path.join(gamma_folder_name, Atoms_folder,atoms_2_name) # This should be the path to the A_{(n+1)/2}B_{(n-1)/2} NP file
    
    # Reading the NPs
    atoms1,atoms2 = read(NP1),read(NP2)
    # Reading the single atom energies
    df_single_atom = pd.read_excel(atom_energies)
    # Reading the NP energies
    df = pd.read_csv(NP_energies)

    def find_base(x):
        """find the base filename from the path

        Args:
            x (str): path to the file

        Returns:
            sol (str): Base filename without the RESTART_ generated by cp2k_helper
        """
        sol = os.path.basename(x)
        if 'RESTART' in sol:
            sol = sol.replace('RESTART_',"")
        return sol

    # Cleaning the data to get the NP names
    df['Folder_Name'] = df['Folder_Name'].apply(find_base) # Get the base folder names

    # Getting the single atom energies
    Pd = float(df_single_atom[df_single_atom['NP']=='Pd']['Energy (Ha)'])*Ha_to_eV
    Pt = float(df_single_atom[df_single_atom['NP']=='Pt']['Energy (Ha)'])*Ha_to_eV
    Au = float(df_single_atom[df_single_atom['NP']=='Au']['Energy (Ha)'])*Ha_to_eV
    Ag = float(df_single_atom[df_single_atom['NP']=='Ag']['Energy (Ha)'])*Ha_to_eV
    Cu = float(df_single_atom[df_single_atom['NP']=='Cu']['Energy (Ha)'])*Ha_to_eV
    d_energies = {'Pd':Pd,'Pt':Pt,'Au':Au,'Ag':Ag,'Cu':Cu}

    # Calculating the NP cohesive energies from DFT (PBE+D3) Calculations
    try: 
        CE_AB1 = ((float(df[df['Folder_Name']==atoms_1_comp]['Energy (eV)']) - (73*d_energies[Atom_Type_1]+74*d_energies[Atom_Type_2]))/147)
        CE_AB2 = ((float(df[df['Folder_Name']==atoms_2_comp]['Energy (eV)']) - (74*d_energies[Atom_Type_1]+73*d_energies[Atom_Type_2]))/147)
    except: # Naming order issue
        atoms_1_comp = Atom_Type_2 + str(math.floor(number_of_atoms/2)) + Atom_Type_1 + str(math.ceil(number_of_atoms/2))
        atoms_2_comp = Atom_Type_2 + str(math.ceil(number_of_atoms/2)) + Atom_Type_1 + str(math.floor(number_of_atoms/2)) 
        CE_AB1 = ((float(df[df['Folder_Name']==atoms_1_comp]['Energy (eV)']) - (73*d_energies[Atom_Type_2]+74*d_energies[Atom_Type_1]))/147)
        CE_AB2 = ((float(df[df['Folder_Name']==atoms_2_comp]['Energy (eV)']) - (74*d_energies[Atom_Type_2]+73*d_energies[Atom_Type_1]))/147)

    # Solving for the gamma values
    solution = calc_gammas([atoms1,atoms2],[CE_AB1,CE_AB2])

    print(solution)
    #  Update gamma dict 
    with open(gamma_values_path) as f:
        gamma_dict = json.load(f)
        updated_gammas = recursive_update(gamma_dict,solution)

    # Write the updated gamma dict to the file
    with open(gamma_values_path, 'w') as f:
        json.dump(updated_gammas, f)