import pickle
import numpy as np
from pathlib import Path
from dolfin import HDF5File, XDMFFile, FunctionSpace, Function, Mesh

class TimeSeriesStorage:
    """ TODO:
            * Add a name-tag for functions 
            * Initialize from function, and get mesh and functionspace from there
            * Consider a switch for turning off xdmf-storage if not needed.
    """
    def __init__(self, mode: str, filepath, *, mesh=None, V=None, name=""):
        self.filepath = Path(filepath).resolve()
        self.name = name
        self.filename = str(self.filepath / self.name)
        if mode == "w":
            assert (mesh is not None and V is not None), 'mode "w" requires a mesh and functionspace V.'
            self._init_write(mesh, V)
        elif mode == "r":
            self._init_read()
        else:
            raise ValueError(f"Invalid mode '{mode}'. Shold be 'w' (write) or 'r' (read).")
    
    def _init_write(self, mesh, V):
        self.mesh = mesh
        self.V = V
        
        # Create directory, value-file and xdmf file for paraview visualization.
        self.filepath.mkdir(exist_ok=True)
        self.hdf5 = HDF5File(mesh.mpi_comm(), str(self.filepath / self.name / f"values.hdf5"), "w")
        self.xdmf = XDMFFile(mesh.mpi_comm(), str(self.filepath / self.name / f"visual.xdmf"))
        
        # Store the mesh to file for simple reconstruction of functionspace. 
        with HDF5File(mesh.mpi_comm(), str(self.filepath / "mesh.hdf5"), "w") as f:
            f.write(mesh, "/mesh")
        
        # Initialize time and length
        self.idx = 0
        self.times = np.empty(0)
        
    def _init_read(self):
        assert self.filepath.is_dir(), f"filepath {self.filepath} does not exist."
        
        # Load mesh from stored file.
        self.mesh = Mesh()
        with HDF5File(self.mesh.mpi_comm(), str(self.filepath / "mesh.hdf5"), "r") as f:
            f.read(self.mesh, "/mesh", True)
            
        # Open file to read values
        self.hdf5 = HDF5File(self.mesh.mpi_comm(), str(self.filepath / self.name / f"values.hdf5"), "r")
        
        # Reconstruct functionspace, and create function to load values.
        self.times, el = self.load_info()
        
        self.V = FunctionSpace(self.mesh, el)
        self.u = Function(self.V)
        
        # Find the number of elements written to file, to determine length.
        self.idx = self.times.size
        

    def _extend_time_vector(self):
        self.times = np.concatenate((self.times, np.nan * np.empty(100)))
            
    def write(self, u, t):
        # Attempt to preallocate data
        try:
            self.times[self.idx] = t
        except IndexError:
            self._extend_time_vector()
            self.times[self.idx] = t
        # Store data
        self.hdf5.write(u.vector(), f"/values_{self.idx}")
        self.xdmf.write(u, t)
        self.idx += 1
        
    def store_info(self):
        with open(self.filepath / "info.pkl", "wb") as f:
            pickle.dump((self.times[~np.isnan(self.times)], self.V.ufl_element()), f)
            
    def load_info(self):
        with open(self.filepath / "info.pkl", "rb") as f:
            times, el = pickle.load(f)
        return times, el

    def close(self):
        
        if hasattr(self, 'hdf5'):
            self.hdf5.close()
        if hasattr(self, 'xdmf'):
            self.xdmf.close()
        
        
    def __iter__(self):
        for i in range(len(self)):
            # u = Function(self.V)

            self.hdf5.read(self.u.vector(), f"/values_{i}", True)
            yield self.u

    def dual_iter(self):
        for i in range(len(self)):
            self.hdf5.read(self.u.vector(), f"/values_{i}", True)
            yield self.times[i], self.u

            
    def __len__(self):
        return self.idx
    
    def __del__(self):
        self.close()
