from ce_expansion.atomgraph.bcm import BCModel
from ce_expansion.atomgraph.adjacency import build_bonds_arr
from ce_expansion.ga.ga import GA
from ce_expansion.ga.ga import Nanoparticle as NP_GA

import ase.cluster as ac
from ase.io import read
from ase.visualize import view

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import sys
import os
from collections import defaultdict
import collections.abc 
import argparse
import json
from IPython.display import HTML,display, Image
import molgif


gamma_folder_name = os.path.join(os.path.dirname(__file__), "Data")
gamma_values_path = os.path.join(gamma_folder_name, "np_gammas.json")

with open(gamma_values_path) as f:
    gammas_np = json.load(f)

ce_bulk_pbe_d3 = {'Au':-3.64,'Pd':-4.20,'Pt':-6.20,"Ag":-2.96,"Cu":-3.95,"Ni":-5.11,"Ir":-7.95} # eV/atom PBE-D3 SOURCE: https://aip.scitation.org/doi/suppl/10.1063/1.4948636/suppl_file/supplementary_material.pdf

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

def get_ordering(atoms):
    """Get the ordering of the atoms in the xyz file

    Args:
        atoms (ase.Atoms): atoms object

    Returns:
        ordering (list): list of the ordering of the atoms in the xyz file
    """

    atom_types = list(np.unique(atoms.symbols))
    atom_types.sort()
    order = []
    for sym in atoms.symbols:
        order.append(atom_types.index(sym))
    return np.array(order)

def get_cutoffs(atoms,x):
    """Custom Cutoffs from custom Radii
        Please add your own custom radii here (covalent radii may be a good place to start so I have implemented that as a default)"""
    radii = {'Au':1.47*x,
             'Pd':1.38*x,
             'Pt':1.38*x}

    # If we have not implemented the radii for the metal in the atoms object, we will use the covalent radii
    df = pd.read_html('http://crystalmaker.com/support/tutorials/atomic-radii/index.html',header=0)[0]
    for element in np.unique(atoms.symbols):
        if element not in radii.keys():
            try:
                radii[element] = float(df[df['ElementSymbol']==element]['CovalentRadius [Å]'].values[0])*x
            except:
                radii[element] = float(df[df['Element Symbol']==element]['Covalent Radius [Å]'].values[0])*x
    return [radii[atom_type] for atom_type in atoms.symbols]

def make_bcm(atoms,x=1.200,CN_Method = 'frac'):
    """Make a BCModel object.  The BCModel object is helpful for calculating the CE of the atoms object as well as to calculate the coordination numbers of the atoms object and finding the shell numbers.  

    Args:
        atoms (ase.Atoms): atoms object
        x (float, optional): scaling factor for the cutoffs. Defaults to 1.200.
        CN_Method (str, optional): Method for calculating coordination number. Defaults to 'frac'.

    Returns:
        bcm (BCModel): BCModel object
    """
    radii = get_cutoffs(atoms,x)
    bonds = build_bonds_arr(atoms,radii)
    bcm = BCModel(atoms,bond_list=bonds,CN_Method=CN_Method)
    # Updating the gamma dictionary with the new gamma values (if new gamma values are available)
    old_gammas = bcm.gammas
    new_gammas = recursive_update(old_gammas,gammas_np)
    bcm.gammas = new_gammas
    ce_bulk_old = bcm.ce_bulk
    ce_bulk_new = recursive_update(ce_bulk_old,ce_bulk_pbe_d3)
    bcm.ce_bulk = ce_bulk_new
    bcm._get_precomps()
    return bcm

def get_comps(atoms,unique_metals):
    """Get the composition of the atoms object

    Args:
        atoms (ase.Atoms): atoms object
        unique_metals (list): list of unique metals in the atoms object
    
    Returns:
        comps (dict): composition of the atoms object
    """
    COMPS = []
    for j in unique_metals:
        COMPS.append(sum(atoms.symbols==j))
    return COMPS

class Nanoparticle:
    def __init__(self,structure,x=1.20,describe="none",method='frac',spike=False):
        """Initialize the Nanoparticle object.  This is a wrapper for the BCModel object and the GA object.  The BCModel object is helpful for calculating the CE of the atoms object as well as to calculate the coordination numbers of the atoms object and finding the shell numbers.  The GA object is helpful for finding the optimal chemical ordering of the atoms object using the BCModel.

        Args:
            structure (str): path to the xyz file or the atoms object
            x (float, optional): scaling factor for the cutoffs. Defaults to 1.200.
            describe (str, optional): description of the nanoparticle. Defaults to "none".
            method (str, optional): Method for calculating coordination number. Defaults to 'frac'.
            spike (bool, optional): Whether or not to spike the GA initial generation with the current NP ordering. Defaults to False.
        """
        # If the xyz file is a string, then read the file, if it is an atoms object, then just use it
        if isinstance(structure,str):
            self.atoms = read(structure)
        else:
            self.atoms = structure
        
        self.cn_method = method # Coordination number method
        self.x = x # Scaling factor for the cutoffs
        self.describe = describe # Description of the nanoparticle
        self.unique_metals = list(np.unique(self.atoms.symbols))
        self.unique_metals.sort()
        self.composition = get_comps(self.atoms,self.unique_metals)
        self.bcm = make_bcm(self.atoms,x=x,CN_Method=method)
        self.bcm_int = BCModel(self.atoms,CN_Method='int')
        self.atom_cut = self.x_cut(self.atoms)
        self.shells,self.comps,self.totals = self.core_shell_info()
        self.df_colors = pd.read_html('https://sciencenotes.org/molecule-atom-colors-cpk-colors/',header=0)[1] # Web-scraping the CPK colors for the atoms
            
        self.GA_init = self.Generate_GA(self.bcm,self.composition,x=x,describe=describe,method=method)
        if spike:
            self.NP_spike = NP_GA(self.bcm,self.composition,get_ordering(self.atoms))
            self.GA_init.pop[0] = self.NP_spike
            self.GA_init.sort_pop()
        
    
    def __len__(self):
        return len(self.atoms)
    
    def core_shell_info(self):
        """Collecting core/shell information from the xyz file

        Returns:
            shells (list): list of shell numbers
            comps (dict): dictionary of compositions for each shell
            totals (list): list of total number of atoms in each shell
        """
        shells = []
        comp = []
        totals = []
        comps = defaultdict(list)
        bcm_int = self.bcm_int
        for i in range(len(bcm_int.shell_map)):
            if i == 0:
                continue
            if i == 1:
                shell_map_core = np.append(bcm_int.shell_map[0],bcm_int.shell_map[i])
                total = len(shell_map_core) 
                totals.append(total)
                for metal_type in self.unique_metals:
                    comps[metal_type].append(sum(self.atoms[shell_map_core].symbols==metal_type)/total)
                shells.append(i)
            else:
                total = len(bcm_int.shell_map[i]) 
                totals.append(total)
                for metal_type in self.unique_metals:
                    comps[metal_type].append(sum(self.atoms[bcm_int.shell_map[i]].symbols==metal_type)/total)
                shells.append(i)
        
        return  shells,comps,totals

    def Generate_GA(self,bcm,COMPS,x=1.20,describe="none",method='frac'):
        return GA(bcm,COMPS,describe)
    
    def x_cut(self,original_atoms):
        atoms = original_atoms.copy()
        core_atom = atoms[self.bcm_int.shell_map[0]][0]
        cutoff = core_atom.a # x coordinate of atom in the core
        atoms_to_del = np.where(atoms.get_positions()[:,0]>cutoff)[0] # Finding all the atom idxs that have a x coord greater than the centers
        del atoms[atoms_to_del]
        return atoms

    def run_ga(self,max_gens=-1,max_nochange=2000):
        """Run the GA to find the optimal chemical ordering.  This function will run the GA until the max number of generations is reached 
            or the max number of generations without a change in the best fitness is reached.

        Args:
            COMPS (dict): dictionary of compositions
            x (float, optional): scaling factor for the cutoffs. Defaults to 1.20.
            describe (str, optional): description for the GA. Defaults to "none".
            method (str, optional): Method for calculating coordination number. Defaults to 'frac'.

        Returns:
            ga (GA): GA object
        """
        ga = self.GA_init
        ga.run(max_gens=max_gens,max_nochange=max_nochange)
        print("Saving optimized structure...")
        self.ga = ga
        self.atoms = self.ga.make_atoms_object()
        self.bcm = make_bcm(self.atoms,x=self.x,CN_Method=self.cn_method)
        self.shells,self.comps,self.totals = self.core_shell_info()
        print("Done!")


    def view(self,cut=False,rotate=False,path=None):
        """View the atoms object

        Args:
            cut (bool, optional): Whether to slice the atoms object.  
            If true a slice of the atoms object in the x-direction will be made before visualizing the structure. 
            Defaults to False.
            rotate (bool, optional): Whether to rotate the atoms object with molgif.
            If true a gif will be created with molgif before visualizing the structure.
            Defaults to False.
            path (str, optional): Path to save the gif. Defaults to None.


        Returns:
            view (ASE): ASE view object
        """
        if cut and not rotate:
            view(self.atom_cut)
        elif cut and rotate:
            if path is None:
                path = "atoms_gif.gif"
            
            if os.path.exists(path):
                os.remove(path)
            molgif.rot_gif(self.atom_cut,optimize=True,save_path=path,overwrite=True,draw_bonds=False,draw_legend=True);
            plt.clf()
            plt.close()
            display(Image(filename=path))
        elif rotate:
            if path is None:
                path = "atoms_full_gif.gif"
                
            if os.path.exists(path):
                os.remove(path)
            molgif.rot_gif(self.atoms,optimize=True,save_path=path,overwrite=True,draw_bonds=False,draw_legend=True);
            plt.clf()
            plt.close()
            display(Image(filename=path))
        else:
            view(self.atoms)
        
    def write(self,filename):
        """Write the atoms object to a file

        Args:
            filename (str): path to the file
        """
        self.atoms.write(filename)

    def calc_ce(self):
        """Calculate the cohesive energy of the nanoparticle

        Returns:
            ce (float): cohesive energy of the nanoparticle
        """
        ce = self.bcm.calc_ce(get_ordering(self.atoms))
        return ce
    
    
    
    def core_shell_plot(self,save=False,saveas='NP_Comp',dpi=300):
        """Plotting core/shell composition on a bar plot
        
        Args:
            save (bool, optional): Whether to save the plot. Defaults to False.
            saveas (str, optional): Name of the file to save the plot as. Defaults to 'NP_Comp'.
            dpi (int, optional): DPI of the saved plot. Defaults to 300.
            
        Returns:
            fig (matplotlib): matplotlib figure object
        """
        
        f, ax = plt.subplots(1, 1, sharey=False)
        ax2 = ax.twiny()
        bottom = np.zeros(len(self.comps[list(self.comps.keys())[0]]))
        for i, key in enumerate(list(self.comps.keys())):
            value = self.comps[key]
            element_color = "#" + self.df_colors[self.df_colors['Element']==key]['Hexadecimal Web Color'].values[0]
            ax.bar(self.shells,value,width=1,bottom=bottom,color= element_color,linewidth=1,edgecolor='black')
            bottom += np.array(value)
        
        surface_num =self.shells[-1]
        ax.set_ylim([0,1])
        ax.set_xlim([0.5,self.shells[-1]+0.5])
        ax.set_xticks(self.shells)
        ax.set_xticklabels(self.shells)
        plt.tight_layout()
        ax.set_xlabel(f'Shell Number (1=Core, {self.shells[-1]}=Surface)',size=20,weight='bold')
        ax.set_ylabel(f'Shell Compositions',size=20,weight='bold')
        ax.legend(self.unique_metals,loc='lower left',fontsize=18)
        ax2.set_xlim(ax.get_xlim()) # ensure the independant x-axes now span the same range
        ax2.set_xticks(self.shells) # copy over the locations of the x-ticks from the first axes
        ax2.set_xticklabels(self.totals) # But give them a different meaning
        # change the fontsize of the xticks and yticks
        ax.tick_params(labelsize=18)
        ax2.tick_params(labelsize=18)
        ax2.set_xlabel(f'Number of Atoms',size=20,weight='bold')
        if save:
            plt.savefig(f'{saveas}.png',dpi=dpi,bbox_inches='tight')
        plt.show()