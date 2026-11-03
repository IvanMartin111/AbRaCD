#!/usr/bin/python3

# iPython interpreter for debugging and testing: 
# /home/mon/prgs/PyRosetta.Ubuntu-12.04LTS.64Bit.monolith/ipython.py

import sys
if '-h' in sys.argv or '--help' in sys.argv or '--help' in sys.argv:
	print('''
This mini-program loads a PDB file as a PyRosetta's pose.

Arguments:
  <pdb_in> <pdb_out> [params]

Optionally, you can provide Rosetta's \".parms\" files for the ligands of the input PDB. 
These files must be provided as a comma-separated list of filenames. 
These kind of files are like those present in directory: 
<Your_PyRosetta>/database/chemical/residue_type_sets/fa_standard/

Example: 
  open_and_close.py 5v71A_Rosetta.pdb 5v71A_Rosetta_curated.pdb GDP.params,8ZG.params

Notes: 
 -PyRosetta will automatically rebuild missing atoms and repack sidechains.
 -Error messages from PyRosetta are very useful to curate and detect PDB file problems.

Mon, 22th September 2017 ''')
	sys.exit(0)

input_pdb = sys.argv[1]
output_pdb = sys.argv[2]
if len(sys.argv) == 4:
	residue_set = sys.argv[3]
	res_set = residue_set.split(',')  # an array is required by generate_nonstandard_residue_set
else:
	residue_set = ""

print("#argv= %d  Input ResidueTypeSet: %s" % (len(sys.argv),residue_set))

# import pyrosetta
#import pyrosetta
#import pyrosetta.rosetta as rosetta
from pyrosetta import init, Pose, pose_from_file

##from pyrosetta.rosetta import *
##rosetta.init()

#import pyrosetta
#import pyrosetta.rosetta as rosetta
#from pyrosetta import init, PyMOLMover, Pose, pose_from_file, ScoreFunction, create_score_function, get_fa_scorefxn, Vector1
#from pyrosetta.rosetta import core, protocols

init()

verb = 'error'
#talaris = " -corrections::restore_talaris_behavior" # required to work with Talaris energies
init(extra_options = "-ex2 -rebuild_disulf false -detect_disulf false -ignore_zero_occupancy false -allow_omega_move false -ignore_unrecognized_res -out:levels core:" + verb + " -out:levels protocols:" + verb )  # + talaris)

print("Loading input pdb: %s" % input_pdb)
pose = Pose() # Create an empty Pose object
if residue_set != "":
	# Read this guide to prepare ligand PDB files for PyRosetta:
	# http://www.pyrosetta.org/obtaining-and-preparing-ligand-pdb-files
	print("Generating a Non-standard residue set")
	restype = generate_nonstandard_residue_set(Vector1( res_set ) )
	# Returns a "custom" ResidueTypeSet with the normal ResidueTypes and any
	#    new ones added as a Vector1 of .params filenames,
	#    the input  <params_list>
	# Example(s):
	#    res_set = generate_nonstandard_residue_set( Vector1( ['ATP.params'] ) )
	print("Using ResidueTypeSet with params: %s" % residue_set)
#	pose_from_pdb(pose, restype, input_pdb) # This automatically rebuilds missing atoms and repacks...
	pose = pose_from_file(restype, input_pdb) # This automatically rebuilds missing atoms and repacks...
else:
#	pose_from_pdb(pose, input_pdb) # This automatically rebuilds missing atoms and repacks...
	pose = pose_from_file(input_pdb) # This automatically rebuilds missing atoms and repacks...


end_res = pose.total_residue()
print("Saving %d residues loaded into output PDB: %s" % (end_res,output_pdb))
pose.dump_pdb(output_pdb)
