#!/usr/bin/env python
import os
import glob
from dataclasses import dataclass
from pathlib import Path

import meshio
import SVMTK as svmtk
from dolfin import Mesh, MeshFunction, XDMFFile, HDF5File, MeshValueCollection
from dolfin.cpp.mesh import MeshFunctionSizet


# TODO: Change all conversions to output the path to resulting files, and create single unified function for conversion.
def prune_z_0(mesh):
    return


def clean_tmp(directory, suffix, no_output=False):
    """ Delete the given directory, with content."""
    dirpath = Path(directory)
    for f in glob.glob(f"{dirpath}/*{suffix}"):
        os.remove(f)

    if len(os.listdir(dirpath)) == 0:
        os.rmdir(dirpath)
    elif no_output:
        return
    else:
        print(f"{dirpath} not empty, and will not be removed.")


def stl2mesh(stlfiles, meshfile_out, resolution, subdomain_map=None, remove_subdomains=None):
    """Creates an svmtk-domain from a set of stl-files, and stores as a meshio .mesh-file. May optionally remove
    sepecifc subdomains, or add subdomain markers."""
    assert type(stlfiles) == list, "stlfiles should be list. (single surface may be wrapped as length 1 list)."
    stlfiles = list(stlfiles)
    surfaces = [svmtk.Surface(str(stl)) for stl in stlfiles]
    if subdomain_map is not None:
        domain = svmtk.Domain(surfaces, subdomain_map)
    else:
        domain = svmtk.Domain(surfaces)
    domain.create_mesh(resolution)
    if remove_subdomains is not None:
        domain.remove_subdomain(remove_subdomains)
    domain.save(str(meshfile_out))


def geo2mesh(infile, outfile, dim=2):
    os.system(f"gmsh -{dim} -format mesh -o {outfile} {infile}")


def meshfunction_default_value(meshfunction, value: int = 0):
    """ Sets the default value for a MeshFunctionSize_t created from a
    MeshValueCollection"""
    for idx, value in enumerate(meshfunction.array() + 1):
        if value == 0:
            meshfunction[idx] = 0
    return meshfunction


def mesh2xdmf(meshfile, xdmfdir, dim=2):
    if dim == 2:
        polytope_label = "triangle"
        facet_label = "line"
    elif dim == 3:
        polytope_label = "tetra"
        facet_label = "triangle"
    else:
        raise ValueError(f"dim should be in (2, 3), got {dim}.")

    # Read meshfile into mesh object, Extract mesh, boundaries and labels.
    mesh = meshio.read(meshfile)
    if dim == 2:
        mesh.prune_z_0()
    points = mesh.points
    polytopes = {polytope_label: mesh.cells_dict[polytope_label]}
    facets = {facet_label: mesh.cells_dict[facet_label]}
    subdomains = {"subdomains": [
        mesh.cell_data_dict["medit:ref"][polytope_label]]}
    boundaries = {"boundaries": [
        mesh.cell_data_dict["medit:ref"][facet_label]]}

    # Write the mesh into new xdmf file
    meshdata = meshio.Mesh(points, polytopes)
    meshio.write(f"{xdmfdir}/mesh.xdmf", meshdata)

    # Write the subdomainds of the mesh
    subdomainfile = meshio.Mesh(points, polytopes, cell_data=subdomains)
    meshio.write(f"{xdmfdir}/subdomains.xdmf", subdomainfile)

    # Write the boundaries/interfaces of the mesh
    boundaryfile = meshio.Mesh(points, facets, cell_data=boundaries)
    meshio.write(f"{xdmfdir}/boundaries.xdmf", boundaryfile)


def xdmf2hdf(xdmfdir, hdf5file):
    # Read xdmf-file into a FEniCS mesh
    dirpath = Path(xdmfdir)
    mesh = Mesh()
    with XDMFFile(str(dirpath / "mesh.xdmf")) as meshfile:
        meshfile.read(mesh)

    # Read cell data to a MeshFunction (of dim n)
    n = mesh.topology().dim()
    subdomains = MeshFunction("size_t", mesh, n)
    with XDMFFile(str(dirpath / "subdomains.xdmf")) as subdomainfile:
        subdomainfile.read(subdomains, "subdomains")

    # Read face data into a Meshfunction of dim n-1
    bdrycollection = MeshValueCollection("size_t", mesh, n-1)
    with XDMFFile(str(dirpath / "boundaries.xdmf")) as boundaryfile:
        boundaryfile.read(bdrycollection, "boundaries")
    boundaries = MeshFunction("size_t", mesh, bdrycollection)
    meshfunction_default_value(boundaries, 0)

    # Write all files into a single h5-file.
    with HDF5File(mesh.mpi_comm(), str(hdf5file), "w") as f:
        f.write(mesh, "/mesh")
        f.write(subdomains, "/subdomains")
        f.write(boundaries, "/boundaries")


def geo2hdf(infile, outfile, dim=2, tmpdir="", clean=True):
    """ Single file for creating h5-file from gmsh .geo file."""
    # Create tmpdir
    inpath = Path(infile)
    tmppath = inpath.parent / tmpdir
    tmppath.mkdir(exist_ok=True)
    meshfile = Path(f"{tmppath}/{inpath.stem}.mesh").resolve()

    # Go from geo->mesh->xdmf->h5
    geo2mesh(inpath, meshfile, dim)
    mesh2xdmf(meshfile, tmppath, dim)
    xdmf2hdf(tmppath, outfile)

    # And clean temporary files created.
    if tmpdir == "/tmp/ratbrain":
        clean_tmp(tmppath, "mesh", no_output=True)
        clean_tmp(tmppath, "xdmf", no_output=True)
        clean_tmp(tmppath, "h5")
    return outfile


def stl2hdf(stlfiles, outfile, resolution, subdomain_map=None, remove_subdomains=None, tmpdir="/tmp/ratbrain"):
    tmppath = create_tmpdir(tmpdir)
    if tmppath != "/tmp/ratbrain":
        print("Intermediate files stored in: ", tmppath)

    stl2mesh(stlfiles, tmppath / 'meshfile.mesh', resolution, subdomain_map=subdomain_map,
             remove_subdomains=remove_subdomains)
    mesh2hdf(tmppath / 'meshfile.mesh', outfile, dim=3, tmpdir=tmpdir)  # Dont think stls will be relevant in 2-dim case

    if tmpdir == "/tmp/ratbrain":
        clean_tmp(tmppath, "meshfile.mesh")


def create_tmpdir(tmpdir):
    tmppath = Path(tmpdir).resolve()
    tmppath.mkdir(exist_ok=True)
    return tmppath


def mesh2hdf(infile, outfile, dim, tmpdir="/tmp/ratbrain"):
    """Single function for creating an h5-file from a meshfile of given dimension.
    TODO: Write test function"""
    # Create tmpdir
    inpath = Path(infile).resolve()
    tmppath = create_tmpdir(tmpdir)
    print(tmppath)
    # meshfile = Path(f"{tmppath}/{inpath.stem}.mesh").resolve()

    # Go from mesh->xdmf->h5
    mesh2xdmf(inpath, tmppath, dim)
    xdmf2hdf(tmppath, outfile)

    # And clean temporary files created.
    if tmpdir == "/tmp/ratbrain":
        clean_tmp(tmppath, "xdmf", no_output=True)
        clean_tmp(tmppath, "h5", no_output=True)
    return outfile


def hdf2fenics(hdf5file, pack=False):
    """ Function to read h5-file with annotated mesh, subdomains
    and boundaries into fenics mesh """
    mesh = Mesh()
    with HDF5File(mesh.mpi_comm(), str(hdf5file), "r") as hdf:
        hdf.read(mesh, "/mesh", False)
        n = mesh.topology().dim()
        subdomains = MeshFunction("size_t", mesh, n)
        hdf.read(subdomains, "/subdomains")
        boundaries = MeshFunction("size_t", mesh, n-1)
        hdf.read(boundaries, "/boundaries")

    if pack:
        return Domain(mesh, subdomains, boundaries)

    return mesh, subdomains, boundaries


@dataclass
class Domain:
    mesh: Mesh
    subdomains: MeshFunctionSizet
    boundaries: MeshFunctionSizet


def unpack_domain(domain: Domain):
    return domain.mesh, domain.subdomains, domain.boundaries
