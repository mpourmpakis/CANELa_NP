import os 
from CANELa_NP.Nanotools import Nanoparticle
import ase.cluster as ac

def test_bimetallic_np_ce():
    atoms = ac.Icosahedron('Au', 5)
    atoms.symbols[100:] = 'Pd'
    NP = Nanoparticle(atoms)
    cohesive_energy_AuPd = round(NP.calc_ce(),12)
    assert cohesive_energy_AuPd == -3.397886376398

def test_trimetallic_np_ce():
    AuPdPt_Path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'Example_Data', 'AuPdPt.xyz'))
    NP = Nanoparticle(AuPdPt_Path)
    cohesive_energy_AuPdPt = round(NP.calc_ce(),12)
    assert cohesive_energy_AuPdPt == -4.0362403098