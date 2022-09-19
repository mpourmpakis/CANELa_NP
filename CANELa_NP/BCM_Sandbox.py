import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from ce_expansion.atomgraph.bcm import BCModel
from ce_expansion.atomgraph.adjacency import build_bonds_arr
from scipy.stats.mstats import gmean,hmean
from ase.io import read
import ase
import ase.neighborlist
from scipy.stats.mstats import gmean,hmean

def get_ordering(atoms):
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

    df = pd.read_html('http://crystalmaker.com/support/tutorials/atomic-radii/index.html',header=0)[0]
    for element in np.unique(atoms.symbols):
        if element not in radii.keys():
            radii[element] = float(df[df['ElementSymbol']==element]['"Crystal"Radius [Å]'].values[0])*x
    return [radii[atom_type] for atom_type in atoms.symbols]

def make_bcm(atoms,x=1.200,CN_Method = 'int'):
    radii = get_cutoffs(atoms,x)
    bonds = build_bonds_arr(atoms,radii)
    bcm = BCModel(atoms,bond_list=bonds,CN_Method=CN_Method)
    return bcm


class BCM_Mod:
    """Original BCM Implementation"""
    def __init__(self,atoms,gammas=False,ce_bulk=False,x=1.200):
        self.atoms = atoms # ASE Atoms object

        self.bcm = make_bcm(atoms,x) # Build the BCM from original code
        
        if not gammas:
            self.gammas = self.bcm.gammas # BCM gammas
        else:
            self.gammas = gammas # Custom gammas
        if not ce_bulk:
            self.ce_bulk = self.bcm.ce_bulk # BCM bulk CE
        else:
            self.ce_bulk = ce_bulk # Custom bulk CE

        self.cns = self.bcm.cn  # coordination numbers
        self.Cb = 12 # Bulk CN - Assuming this is a FCC metal

        self.syms = atoms.symbols # atom symbols
        self.atom_types = np.sort(np.unique(atoms.symbols)) # unique atom types


    def calc_ce(self):
        """Calculate the CE of the metal nanoparticle with the modified BCM"""
        num_sum = 0
        for i,j in self.bcm.bond_list:
            A = self.atoms.symbols[i]
            B = self.atoms.symbols[j]
            part_1 = self.gammas[A][B]*(self.ce_bulk[A]/self.cns[i])*np.sqrt(self.cns[i]/self.Cb) 
            part_2 = self.gammas[B][A]*(self.ce_bulk[B]/self.cns[j])*np.sqrt(self.cns[j]/self.Cb) 
            num_sum += (part_1 + part_2)
        return (num_sum/(len(self.atoms)*2)) # CE of the metal nanoparticle
    