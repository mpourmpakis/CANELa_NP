<h1 align="center">CANELa_NP</h1>

Python Package developed for Demystifying the Chemical Ordering of Multimetallic Nanoparticles by Dennis Loevlie, Brenno Ferreira, and Giannis Mpourmpakis*

-------
# Installation 

```
git clone https://github.com/mpourmpakis/CANELa_NP.git
cd CANELA_NP
pip install -e .
```

<h1 align="center">Package Tutorial (demonstration of use)</h1>


## Importing packages

```python
from CANELa_NP.Nanotools import Nanoparticle
import ase.cluster as ac
```

## Creating a bimetallic NP with ASE 


```python
atoms = ac.Icosahedron('Au',5) 
```


```python
print(atoms)
```




    >>> Atoms(symbols='Au309', pbc=False, tags=...)




```python
atoms.symbols[100:] = 'Pd'
```


```python
print(atoms)
```




    >>> Atoms(symbols='Au100Pd209', pbc=False, tags=...)



## Initializing the NP Object


```python
NP = Nanoparticle(atoms)
```

# Visualizing the core to shell element distribution


```python
NP.core_shell_plot()
```


    
![png](README_Notebook_files/README_Notebook_10_0.png)
    


## Optimizing the chemical ordering


```python
NP.run_ga(max_gens=-1,max_nochange=1_000)
```

    --------------------------------------------------
    GA Sim for Au100Pd209 - none:
     Min: -3.66177 eV/atom -- Gen: 02840
     Form: Au100Pd209
    nAtom: 309
    nGens: 2840
    Start: -3.44202 eV/atom
     Best: -3.66177 eV/atom
     Diff: -0.21974 eV/atom (6.38%)
     Mute: 80.0%
      Pop: 50
     Time: 0:00:28
    --------------------------------------------------
    Saving optimized structure...
    Done!


## Visualizing the core/shell distribution of the optimized chemical ordering


```python
NP.core_shell_plot()
```


    
![png](README_Notebook_files/README_Notebook_14_0.png)
    


## Visualizing the NP with ase gui

Viewing the full NP


```python
NP.view()
```

![png](README_Notebook_files/full_np.png)


Viewing a slice of the NP (in the x-direction)


```python
NP.view(cut=True)
```

![png](README_Notebook_files/half_np.png)

## Working with your own xyz files


```python
xyz_file = 'Example_Data/AuPdPt.xyz'
```


```python
NP = Nanoparticle(xyz_file)
```


```python
NP.core_shell_plot()
```


    
![png](README_Notebook_files/README_Notebook_23_0.png)
    

