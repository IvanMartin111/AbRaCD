#!/usr/bin/env python



####################################################################################
# MAIN PYTHON CODE
####################################################################################

# PARSER
# ----------------------------------------------------------------------------------
# parser object for managing input options
# all defaults are for the example using "test_in.pdb"
import optparse    # for option sorting
usage = "USAGE: %prog [options]"
parser = optparse.OptionParser(usage=usage)
parser.add_option('--base_in', dest = 'base_in',
    default = '',    # default to False
    help = 'Base-name of input files: <base>.pdb , <base>_dh.txt , <base>_val.txt , <base>_len.txt' )
parser.add_option('--pdb_in', dest = 'pdb_in',
    default = '',    # default example PDB
    help = 'File name of the PDB containing the minimal backbone coordinates (N,CA,C) of the loop to remodel' )
parser.add_option('--base_out', dest = 'base_out',
    default = 'rcd',    # default to False
    help = 'Output base-name (default=rcd)' )
parser.add_option('--ri', dest = 'ri',
    default = '',    # default to False
    help = 'Nt-anchor residue index of the loop (PDB numbering)' )
parser.add_option('--rf', dest = 'rf',
    default = '',    # default to False
    help = 'Ct-anchor residue index of the loop (PDB numbering)' )
parser.add_option('--chain', dest = 'chain',
    default = '',    # default to False
    help = 'Chain ID (default = chain ID of the first residue)' )
parser.add_option('--Bfactor', dest = 'Bfactor',
    default = 50,    # default scale factor for B-factors
    help = 'Scale factor for B-factors (default=50)' )
parser.add_option('--nmin', dest = 'nmin',
    default = 0,    # default to False
    help = 'Number of minimization steps of the Monte-Carlo plus Minimization algorithm of Li & Scheraga (default=0)' )
parser.add_option('--normsds', dest = 'normsds',
    default = '',    # default None
    help = 'Disables Rosetta\'s RMSDs computation (default = disabled)' )
parser.add_option('--initneighbor', dest = 'initneighbor',
    default = 0,    # default to False
    help = '>0 --> CB-CB distance to define loop neighborhood for initial repacking and minimization, 0 --> loop-only, <0 --> no-repacking/minimization (default=-1).' )
parser.add_option('--nminneigh', dest = 'nminneigh',
    default = '',    # default to False
    help = 'Number of loops used to define --initneighbor neighborhood. Warning: <base>_dh.txt <base>_val.txt and <base>_len.txt files required. By default, the loop already present in PDB will be used. (default=disabled)' )
parser.add_option('--repackouther', dest="repackouther", 
    default=False, # default inner loop
    action="store_true",
    help = 'Repack loop in outher loop (in the \"--irepack\" selected iteration) (default=disabled, inner loop repacking)' )
parser.add_option('--irepack', dest = 'irepack',
    default = -1,    # default, disabled
    help = 'Inner or Outher loop iteration to trigger loop-repacking. -1= disabled, 0= all, 1= first, etc... (default=-1, disabled)' )
parser.add_option('--mintype', dest = 'mintype',
    default = 'dfpmin_armijo_nonmonotone',    # default to False
    help = 'Flavors of descent minimization in PyRosetta. You can select among any \"Min mover\" minimization method of Rosetta: dfpmin (exact line search), dfpmin_armijo (inexact line search), dfpmin_armijo_nonmonotone (less exact line search), dfpmin_strong_wolfe (More-Thuente line minimization). See: https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d7/df1/minimization_overview.html for reference.' )
parser.add_option('--energy', dest = 'energy',
    default = 'talaris2014',    # default energy function
    help = 'Rosetta energy function: talaris2013, talaris2014, etc...' )
parser.add_option('--fas', dest = 'fas',
    default = '',    # default to False
    help = 'VdW repulsive factors for fast-relax minimization protocol (default= 0.05,0.30,0.60,1.00)' )
parser.add_option('--tols', dest = 'tols',
    default = '',    # default to False
    help = 'Convergence/tolerance factors for quick-relax minimization protocol (default= 0.01,0.01,0.001,0.00001)' )
parser.add_option('--ntpsi', dest="ntpsi", 
    default=False, 
    action="store_true",
    help = 'Enable N-terminal-PSI dihedral angle optimization. (default=disabled)' )
parser.add_option('--deleteloop', dest="deleteloop", 
    default=False, 
    action="store_true",
    help = 'Delete loop (without anchors) before anything. (default=disabled)' )
parser.add_option('--cure', dest="cure", 
    default=False, 
    action="store_true",
    help = 'PyRosetta\'s automatic curation. Just load with pose_from_pdb() and write PDB file. (default=disabled)' )

(options,args) = parser.parse_args()

# os.environ['LD_LIBRARY_PATH']='/home/mon/prgs/PyRosetta.Ubuntu-12.04LTS.64Bit.monolith/'

import time
start = time.time()
import os;

import pyrosetta
import pyrosetta.rosetta as rosetta
from pyrosetta import init, Pose, pose_from_file, ScoreFunction, create_score_function #, get_fa_scorefxn, Vector1



# Import my stuff into current namespace...
from monfuncs import *

 
# WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION

cure = options.cure

if not cure:
	print('PDB curation disabled')
	init(extra_options = "-ex2 -rebuild_disulf false -detect_disulf false -ignore_zero_occupancy false -allow_omega_move false -ignore_unrecognized_res -corrections::restore_talaris_behavior")
	cure = False
else:
	print('PDB curation enabled (only this task will be performed). Here, -detect_disulf option must be false, otherwise PyRosetta dumps error...')
	init(extra_options = "-ex2 -rebuild_disulf false -detect_disulf false -ignore_zero_occupancy false -allow_omega_move false -ignore_unrecognized_res -corrections::restore_talaris_behavior")
	cure = True

# About repacking: https://www.rosettacommons.org/manuals/rosetta3.1_user_guide/opt_packing.html
# -packing:ex1:level                Use extra chi1 sub-rotamers for all residues that pass the extrachi_cutoff. [Integer]
#                                   (level 1 Default  +/- one standard deviation (sd); 3 samples)
# -packing:extrachi_cutoff          Number of neighbors a residue must have before extra rotamers are used. default: 18 [Integerx]
# -packing:use_input_sc             Use rotamers from input structure in packing. By default, input sidechain coords are NOT
#                                   included in rotamer set but are discarded before the initial pack; with this flag, the input
#                                   rotamers will NOT be discarded. Note that once the starting rotamers are replaced by any mechanism,
#                                   they are no longer included in the rotamer set. (rotamers included by coordinates)
#
# Mon: without omega (like loops_energy7.py) "OMEGA true SEEMS TO BE BAD... so disable it"
# -ex1 option is in "$PYROSETTA/rosetta/__init__.py" file
# See: ~/PROCESS/rcd/loop_modeling/input/preparation/pre_min_pack.py (Rosetta's pre-protocol) (They use -min_pack, check it???)
# Downloaded from: https://guybrush.ucsf.edu/benchmarks/captures/loop_modeling

# -detect_disulf true --> crashes with "1ezm" and "1lst" (8 aa) "ERROR: ! found_aa_difference"
# -detect_disulf true --> Crashes in "task_packer", it delete CYS .... (check another PyRosseta version????) (Fixed: using "vector1" in repacking)

end = time.time()
print("Init time: %f" % (end - start))


# Loading parser variables

base_in = options.base_in

pdb_in = options.pdb_in
if not pdb_in:
	pdb_in = "%s.pdb" % base_in # Take input pdb file name from input base-name

if not pdb_in and not base_in:
	print("ERROR: You must supply an input PDB file or valid Basename!")
	exit()




nminneigh = options.nminneigh
if nminneigh:
	nminneigh = int(nminneigh)
	print("Using <base>_dh.txt <base>_val.txt and <base>_len.txt files to generate loops and compute neighborhood.")
	dh_file = "%s_dh.txt" % base_in # Take input dihedrals file name from input base-name
	val_file = "%s_val.txt" % base_in # Take input valence angles file name from input base-name
	len_file = "%s_len.txt" % base_in # Take input bond lengths file name from input base-name
	# Reading Internal Coordinates
	print("Reading Internal Coordinates from %s, %s, and %s files." % (dh_file,val_file,len_file))
	dihedrals = read_dihedrals(dh_file) # Read dihedral angles from file
	valences = read_valences(val_file) # Read valence angles from file
	lengths = read_lengths(len_file) # Read bond lengths from file
	if nminneigh > int(len(dihedrals)):
		nminneigh = int(len(dihedrals))
		print("Warning: The number of loops to compute neighborhood is higher than the number of loops provided (%d). Setting it to the maximum available." % nminneigh)



base_out = options.base_out
pdb_read = "%s_input.pdb" % base_out # readed pdb with all atoms added by PyRosetta
energy_file = "%s_min.txt" % base_out # Output table data (from Rosetta)

ri = int(options.ri)
if not ri or ri < 0:
	ri = 0;
else:
	ri = int(ri)
print("Parser: ri= %d" % ri)
rf = int(options.rf)
if not rf or rf < 0:
	rf = 0;
else:
	rf = int(rf)
print("Parser: rf= %d" % rf)
chain = options.chain
if chain:
	print("Parser: chain= %c" % chain)

nmin = int(options.nmin)
print("Parser: nmin= %d" % nmin)

initneighbor = options.initneighbor
if not initneighbor:
	initneighbor = int(-1)
else:
	initneighbor = float(initneighbor)
print("Parser: initneighbor= %f" % initneighbor)

repackouther = options.repackouther
if repackouther:
	print("Parser: Loop-Repacking enabled in OUTHER loop.")
else:
	print("Parser: Loop-Repacking enabled in INNER loop.")
irepack = int(options.irepack)
if irepack >= 0:
	if irepack == 0:
		print("Parser: Loop-Repacking enabled always (irepack=%d)" % irepack)
	else:
		print("Parser: Loop-Repacking enabled in iteration %d" % irepack)

else:
	print("Parser: Loop-Repacking disabled (irepack=%d)" %irepack)
	irepack = int(-1); # Repacking disabled

opt_fas = options.fas
if not opt_fas:
	opt_fas = '0.05,0.30,0.60,1.00'
else:
	opt_fas = options.fas
fas=[float(x) for x in opt_fas.split(',')]  # split a comma-separated string into array
print("Using the following VdW factors for quick-relax minimization:")
print(fas)

opt_tols = options.tols
if not opt_tols:
	opt_tols = '0.01,0.01,0.001,0.00001'
else:
	opt_tols = options.tols
tols=[float(x) for x in opt_tols.split(',')]  # split a comma-separated string into array
print("Using the following tolerances for quick-relax minimization:")
print(tols)

mintype = options.mintype
if not mintype:
	mintype = 'dfpmin_armijo_nonmonotone'
print("Using \"%s\" flavor in descent minimization (mintype option for MinMover)." % mintype)

energy = options.energy
if not energy:
	energy = 'talaris2014'
print("Using Rosetta's %s energy function (for repacking and minimization)." % energy)


ntpsi = options.ntpsi
print('Loop N-terminal PSI angle optimization: ', ntpsi)

deleteloop = options.deleteloop
print('Delete loop: ', deleteloop)

# RMSD masks
CA_rmsd = []
CA_rmsd.append('CA')
BB_rmsd = []
BB_rmsd.append('N')
BB_rmsd.append('CA')
BB_rmsd.append('C')
BBO_rmsd = []
BBO_rmsd.append('N')
BBO_rmsd.append('CA')
BBO_rmsd.append('C')
BBO_rmsd.append('O')
BBOCB_rmsd = []
BBOCB_rmsd.append('N')
BBOCB_rmsd.append('CA')
BBOCB_rmsd.append('C')
BBOCB_rmsd.append('O')
BBOCB_rmsd.append('CB')
HA_rmsd = []
HA_rmsd.append('H') # considers only non-Hydrogen atoms
SC_rmsd = []
SC_rmsd.append('H') # non-Hydrogen atoms
SC_rmsd.append('N ') # non N
SC_rmsd.append('CA') # non CA
SC_rmsd.append('C ') # non C
SC_rmsd.append('O ') # non O


#######################################
# MAIN
#######################################

# Load the data from pdb_file into the pose
print("Loading input pdb: %s" % pdb_in)
start = time.time()
pose = Pose() # Create an empty Pose object
# Bug 20/08/2018 in PyRosetta versions above v72...
pose = pose_from_file(pdb_in) # Warning, this automatically rebuilds missing atoms and repacks...
end = time.time()
print("Load time: %f" % (end - start))

# Show pdb info
info = pose.pdb_info()
print(info)

# Write input pdb (now full atom and packed)
pose.dump_pdb(pdb_read)
print("Input PDB file with added missing atoms and hydrogens written: %s\n" % pdb_read)
if cure:
	print("PyRosetta should have added missing atoms and hydrogens to it... -detect_disulf false. Deliberate exit!\n")
	exit()


# Find chain ID if not provided by parser
chain_num = 1
if not chain or len(chain) != 1:
	chain = info.chain(chain_num) # Get chain ID of first residue
	print("Chain %c selected by default" % chain)
else:
	print("Chain %c selected" % chain)

# Find correspondence between PDB and Rosetta numbering schemes for residues
if ri > 0 and rf > 0:
	rri = info.pdb2pose(chain,ri) # initial Rosetta residue index
	print("Initial Rosetta residue index (Nt-anchor) %d" % rri)
	if rri == 0:
		print("Residue %d from chain %c not found! Forcing exit!" % (ri, chain))
		exit()
	rrf = info.pdb2pose(chain,rf) # final Rosetta residue index
	print("Final Rosetta residue index (Ct-anchor) %d" % rrf)
	if rrf == 0:
		print("Residue %d from chain %c not found! Forcing exit!" % (rf, chain))
		exit()

	if (rrf-rri) != (rf-ri):
		print("WARNING: Correspondence mismatch between PDB and Rosetta numbering! (rrf=%d rri=%d rf=%d ri=%d)" % (rrf,rri,rf,ri))

	print("PDB numbered segment %d-%d corresponds to %d-%d in Rosetta." % (ri, rf, rri, rrf))

	if rrf-rri > 1:
		loop_sel = range( rri+1 , rrf ) # loop range selection (for RMSD computations) NOT-INCLUDING ANCHORS
		print("Selected residues range for RMSDs: from %d to %d (loop without anchors)" % (rri+1,rrf-1))  # without anchors?
	else:
		print("WARNING: All loop residues missing!")

else:
	rri = 1;
	rrf = pose.total_residue();
	loop_sel = range( rri , rrf+1 ) # Full protein selection (for RMSD computations) 
	print("Selected residues range for RMSDs: from %d to %d (full protein)" % (rri,rrf))
	print(loop_sel)

# Create a standard ScoreFunction
scorefxn = create_score_function(energy) # state-of-the-art Rosetta's energy

# MON FIXED (12/1/2016)
scorefxn.set_weight(core.scoring.chainbreak,100.0) # 100 (instead of 1) is important to avoid chain breaks
scorefxn.set_weight(core.scoring.coordinate_constraint,1.0)
print(scorefxn)

# Compute native loop Energy (initial pose should have native loop)
poseE0 = scorefxn(pose)

# This trick backups the loop PDBinfo from input molecule
poseEloop = 0.0

dummy = Pose()
dummy.assign(pose) # full copy of "pose"
# Deleting non-loop residues 
# (this seems the only way to assing later a valid PDBinfo object to each loop...)
print("Deleting residues from %d to %d in dummy" % (rrf+1,pose.total_residue()))
rosetta.protocols.grafting.delete_region(dummy,rrf+1,pose.total_residue())
print("Deleting residues from %d to %d in dummy" % (1,rri-1))
rosetta.protocols.grafting.delete_region(dummy,1,rri-1)
loop_info = dummy.pdb_info() # Store valid PDBinfo to be asigned later to each generated loop
poseEloop = scorefxn(dummy) # Get just-loop energy (including anchors)

nat_name = "%s_nat.pdb" % (base_out)
dummy.dump_pdb(nat_name)
print("Loop (with anchors) saved in: %s\n" % (nat_name))

# INITIAL STUFF...
if nminneigh:
	print("NEIGHBORHOOD CALCULATION FROM the first %d SAMPLED LOOPS" % nminneigh)
	start = time.time()
	# Create a full copy the initial pose for every loop (Avoid repacking errors between different loops)
	pose2 = Pose()
	pose2.assign(pose)
	initneighbors = 0  # Neighbors' hash is first initializaed with 0
	for l in range(0,nminneigh):

		# Update Internal Coordinates from loaded ICs arrays
		update_dihedrals(pose2, dihedrals[l], rri)
		update_valences(pose2, valences[l], rri)
		update_lengths(pose2, lengths[l], rri)

		# Neighbors' hash is incrementally built...
		initneighbors = neighbors(pose2,initneighbor,rri,rrf,initneighbors)
		if l-1 % 10 == 0 or l == nminneigh-1:
			print("Neighborhood from loop %d = %d resiudes\n" % (l+1,len(initneighbors)))
	end = time.time()
	print("Initial Neighborhood detection time: %f s" % (end-start))
else:
	if ri > 0 and rf > 0:
		initneighbors = neighbors(pose,initneighbor,rri,rrf)  # Define initial neighborhood for initial pose
	else:
		initneighbors = neighbors(pose)  # Define initial neighborhood for initial pose

if deleteloop:
	# Delete loop residues (not-anchors)
	pose2 = Pose()
	pose2.assign(pose) # Create a full copy the initial pose for RMSDs
	print("Number of residues before deletion: %d " % pose2.total_residue())

	print("Deleting residues from %d to %d in pose before any minimization" % (rri+1,rrf-1))
	rosetta.protocols.grafting.delete_region(pose2,rri+1,rrf-1) # detect_disulfides() seems not neccessary now...

	# This gave Disulphide-bond ERRORS....	Check this forum: https://www.rosettacommons.org/node/3932
	# Like this one: core.conformation.util: ERROR: Couldn't find a disulfide-bonded equivalent for residue CYS28.


	






	print("Number of residues after deletion: %d " % pose2.total_residue())

#	# This is the "third" solution


	# Translating neighborhood selection mask from the initial (complete) protein into the loop-less version.
	print("Current neighborhood in Rosseta numbering scheme")
	print(initneighbors)

	# Translate selection mask from Rosetta indices to PDB indices
	#  neighs --> hash of residue indices in Rosetta numeration scheme
	#  pose --> Pose to get PDB indices from
	# OUTPUT:  hash of residue indices in PDB numeration scheme
	print("Translation from Rosseta numbering scheme to PDB (chain %c)" % (chain))
	pdbneighs = rosseta2pdb(initneighbors,pose)
	print(pdbneighs)

	# Translate selection mask from PDB indices to Rosseta indices
	#  neighs --> hash of residue indices in PDB numeration scheme
	#  pose --> Pose to get Rosseta indices from
	#  chain --> Chain-ID (one-letter) required
	# OUTPUT:  hash of residue indices in Rosetta numeration scheme
	print("Translation from PDB numbering scheme to Rosetta (chain %c)" % (chain))
	initneighbors =	pdb2rosseta(pdbneighs,pose2,chain)
	print(initneighbors)


# Watch out this!!! <-- This seems the problem later...
	loop_sel = initneighbors

else:
	pose2 = pose
	print("Using the following neighborhood for repacking and minimization (Rosseta numbering scheme)")
	print(initneighbors)

# Native Pose backup (to compute RMSDs)
posen = Pose()
posen.assign(pose2) # Create a full copy the initial pose for RMSDs

print("Repacking Initial neighborhood (initneighbor=%f)\n" % (initneighbor))
neighbor_packer(pose2,scorefxn,initneighbors)  # Repack current pose
print("Minimization with Initial neighborhood (initneighbor=%f)\n" % (initneighbor))
if ri > 0 and rf > 0 and not deleteloop:
	monrelax(pose2,scorefxn,nmin,-1,initneighbors,fas,tols,rri,rrf,repackouther,irepack,ntpsi,mintype) # "nmin" iters with all "fas" (slower)
else:
	monrelax(pose2,scorefxn,nmin,-1,initneighbors,fas,tols,0,0,repackouther,irepack,ntpsi,mintype) # "nmin" iters with all "fas" (slower)

# Write minimized pdb
initmin_name = "%s_min.pdb" % (base_out)
pose2.dump_pdb(initmin_name)
print("Minimized PDB structure saved in: %s\n" % (initmin_name))

# Compute Energy of the Minimized initial structure (initial pose should have native loop or not)
poseEmin = scorefxn(pose2)

# Compute RMSDs of the Minimized initial structure
ca_rmsd = rmsd(posen,pose2,loop_sel,CA_rmsd)
bb_rmsd = rmsd(posen,pose2,loop_sel,BB_rmsd)
bbo_rmsd = rmsd(posen,pose2,loop_sel,BBO_rmsd)
bbocb_rmsd = rmsd(posen,pose2,loop_sel,BBOCB_rmsd)
ha_rmsd = rmsd(posen,pose2,loop_sel,None,HA_rmsd)
sc_rmsd = rmsd(posen,pose2,loop_sel,None,SC_rmsd)
full_rmsd = rmsd(posen,pose2,loop_sel,None,None)

print("#%4s %8s %8s %6s %5s %5s %5s %8s %5s %5s %5s %5s %5s %5s %5s" % ('loop', 'E_full', 'E_loop', 'R_rcd', 'R_BB', 'R_HA', 'R_SC', 'E_full2', 'R_CA', 'R_BB', 'R_BBO', 'R_CB', 'R_HA', 'R_SC', 'R_ALL'))
print("#%4s %8.2f %8.2f %6.3f %5.2f %5.2f %5.2f %8.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n" % (' ', poseE0, poseEloop, 0.0, 0.0, 0.0, 0.0, poseEmin, ca_rmsd, bb_rmsd, bbo_rmsd, bbocb_rmsd, ha_rmsd, sc_rmsd, full_rmsd))

f=open(energy_file, 'w+')
f.write("#%4s %8s %8s %6s %5s %5s %5s %8s %5s %5s %5s %5s %5s %5s %5s\n" % ('loop', 'E_full', 'E_loop', 'R_rcd', 'R_BB', 'R_HA', 'R_SC', 'E_full2', 'R_CA', 'R_BB', 'R_BBO', 'R_CB', 'R_HA', 'R_SC', 'R_ALL') )
f.write("#%4s %8.2f %8.2f %6.3f %5.2f %5.2f %5.2f %8.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n" % (' ', poseE0, poseEloop, 0.0, 0.0, 0.0, 0.0, poseEmin, ca_rmsd, bb_rmsd, bbo_rmsd, bbocb_rmsd, ha_rmsd, sc_rmsd, full_rmsd) )


# Write minimized pdb with native loop

exit()


