# MultiCompartment Solute Transport
This repository contains source code for the simulations pressented in the paper "Multi-compartmental model of glymphatic clearance of solutes in brain tissue" by Poulain, Riseth and Vinje. For now it contains only a minimal working example for running the simulations, but will be reorganized to be more accessible and more easily **extended** in the near future. 

## Files
* `main.py`: Run simulations with the multi-compartment model. Parameters, model choices and other simulation settings are all set within the script, and have to be adjusted manually if you want to see changes.
* `ts_storage.py`: Contains a class which wraps some FEniCS storage types. Stores the solutions both as XDMF (for easy visualization in paraview) and HDF5 together with various python objects, for easy loading and post-processing of solutions.

...



## Setup
We recommend using conda to install necessary packages. A working environment may be created using:
```
conda create -n multicomp -c conda-forge fenics matplotlib jupyter
conda activate multicomp
pip install h5py
pip install meshio
pip install git+https://github.com/SVMTK/SVMTK.git
pip install -e . 
```

### Mesh Creation
NB: The information below is not currently correct, but will be adjusted in the future. For now, the meshes are already available as `h5`-files available in the `mesh` directory.
---

Necessary files for recreating the mesh are available within the `mesh/`-folder. The script `mesh_generation.py` reads the `stl`-files within the `mesh/stl-files` directory, and generates the meshes. It may be run by
```
./mesh_generation.py --resolution 32
```

Meshes of varying resolution may be created from a space-separated list of values. The scripts outputs files of format `.h5`, which may be loaded into fenics.

