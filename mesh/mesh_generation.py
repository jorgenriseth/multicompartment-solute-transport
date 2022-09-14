#!/usr/bin/env python
from argparse import ArgumentParser
from pathlib import Path
from meshprocessing import stl2hdf
import SVMTK


parser = ArgumentParser(
    description="""
                --------------
                Convert a pair of stl-files ('brain.stl', 'ventricles.stl') representing the brain and ventricles of
                a ratbrain into a dolfin mesh, where the ventricles are excluded. Assumes the stl-files are stored in
                the same directory, given by the '--inputdir'-argument."""
)
meshdir = f"{Path(__file__).resolve().parent}"
parser.add_argument('--resolution', type=int, required=True, nargs="+")
parser.add_argument('--inputdir', type=str, default=f"{meshdir}/stl-files/")
parser.add_argument('--outputdir', type=str, default=f"{meshdir}")
parser.add_argument('--tmpdir', type=str, default='/tmp/ratbrain/', help='Defaults to "/tmp/ratbrain". If provided,'
                                                                         'any intermediate files required in the '
                                                                         'creation of the mesh file will be saved in '
                                                                         'this directory.')

args = parser.parse_args()

# Parse paths
inputdir = Path(args.inputdir).resolve()
outputdir = Path(args.outputdir).resolve()

# TODO: Allow for differently named output-files? Need to decide how to combine with multiple resolution inputs.
# if args.outfile != meshdir:
#     outfile = Path(args.outfile).resolve()
#     assert outfile.suffix == "h5", "output file should be of type .h5"
#     assert length(args.resolution)

surfaces = [inputdir / "brain.stl", inputdir / "ventricles.stl"]

# Define subdomains
smap = SVMTK.SubdomainMap()
smap.add("10", 1)  # Brain matter
smap.add("01", 2)  # Ventricles
smap.add("11", 2)  # Possible overlap between the two are labelled as ventricles.

for res in args.resolution:
    print()
    print(f"Creating mesh with resolution {res}.")
    output = outputdir / f"mesh{res}.h5"
    stl2hdf(surfaces, output, res, subdomain_map=smap, remove_subdomains=2, tmpdir=args.tmpdir)
