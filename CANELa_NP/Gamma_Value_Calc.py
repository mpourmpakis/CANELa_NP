import pandas as pd
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

# Constants 
"""Please add your own Bulk CE values here for whatever metals you are using"""
ce_bulk_pbe_d3 = {'Au':-3.64,'Pd':-4.20,'Pt':-6.20} # eV/atom
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

    return {A:{A:1,B:float(sol[gamma_a])},B:{B:1,A:float(sol[gamma_b])}}

if __name__ == '__main__':
    gamma_folder_name = os.path.join(os.path.dirname(__file__), "Data")
    gamma_values_path = os.path.join(gamma_folder_name, "np_gammas.json")
    atom_energies = os.path.join(gamma_folder_name, "single_atom_energies.xlsx")
    NP_energies = os.path.join(gamma_folder_name, "NP_energies.csv")
    NP1 = os.path.join(gamma_folder_name, "PtPd","Pd73Pt74.xyz")
    NP2 = os.path.join(gamma_folder_name, "PtPd","Pd74Pt73.xyz")
    atoms1,atoms2 = read(NP1),read(NP2)
    df_single_atom = pd.read_excel(atom_energies)
    df = pd.read_csv(NP_energies)

    def find_base(x):
        sol = os.path.basename(x)
        if 'RESTART' in sol:
            sol = sol.replace('RESTART_',"")
            print(sol)
        return sol

    df['Folder_Name'] = df['Folder_Name'].apply(find_base) # Get the base folder names



    Pd = float(df_single_atom[df_single_atom['NP']=='Pd']['Energy (Ha)'])*Ha_to_eV
    Pt = float(df_single_atom[df_single_atom['NP']=='Pt']['Energy (Ha)'])*Ha_to_eV
    Au = float(df_single_atom[df_single_atom['NP']=='Au']['Energy (Ha)'])*Ha_to_eV

    # print(df[df['Folder_Name']=='Au73Pt74']['Energy (eV)'])
    CE_AB1 = ((float(df[df['Folder_Name']=='Pd73Pt74']['Energy (eV)']) - (73*Pd+74*Pt))/147)
    CE_AB2 = ((float(df[df['Folder_Name']=='Pd74Pt73']['Energy (eV)']) - (74*Pd+73*Pt))/147)
    solution = calc_gammas([atoms1,atoms2],[CE_AB1,CE_AB2])

    
    #  Update gamma dict 
    with open(gamma_values_path) as f:
        gamma_dict = json.load(f)
        updated_gammas = recursive_update(solution,gamma_dict)


    with open(gamma_values_path, 'w') as f:
        json.dump(updated_gammas, f)