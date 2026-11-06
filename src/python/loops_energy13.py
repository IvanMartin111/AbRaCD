#!/usr/bin/env python

##!/usr/bin/python3


####################################################################################
# MAIN PYTHON CODE
####################################################################################
import mpi4py

from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = MPI.COMM_WORLD.Get_size()
print("PyRosetta MPI worker %d awoke from %d!\n" % (rank,size));


# PARSER
# ----------------------------------------------------------------------------------
# parser object for managing input options
# all defaults are for the example using "test_in.pdb"
import optparse    # for option sorting
usage = "usage: %prog [options]"
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
parser.add_option('--dihedrals', dest = 'dh_file',
    default = '',    # default to False
    help = 'Input dihedral angles table in this plain-text format: Omega_1, Phi_1, Psi_1, Omega_2, ..., Omega_N, Phi_N, Psi_N' )
parser.add_option('--valences', dest = 'val_file',
    default = '',    # default to False
    help = 'Input valence angles table in this plain-text format: NCAC_1, CACN_1, CNCA_1, NCAC_2, ..., NCAC_3, CACN_3, CNCA_3, NCAC_N.' )
parser.add_option('--lengths', dest = 'len_file',
    default = '',    # default to False
    help = 'Input bond lengths table in this plain-text format: N-CA_1, CA-C_1, C-N_1, N-CA_2, ..., N-CA_3, CA-C_3, C-N_3, N-CA_4, CA-C_4.' )
parser.add_option('--rmsds', dest = 'rmsds_file',
    default = '',    # default to False
    help = 'Min-RMSDs file name' )
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
parser.add_option('--loopI', dest = 'loopI',
    default = '',    # default None
    help = 'Initial loop to be procesed (default = 1' )
parser.add_option('--loopF', dest = 'loopF',
    default = '',    # default None
    help = 'Final loop to be procesed (default = last loop index)' )
parser.add_option('--normsds', dest = 'normsds',
    default = '',    # default None
    help = 'Disables Rosetta\'s RMSDs computation (default = disabled)' )
parser.add_option('--initneighbor', dest = 'initneighbor',
    default = 0,    # default to False
    help = '>0 --> CB-CB distance to define loop neighborhood for initial repacking and minimization, 0 --> loop-only, <0 --> no-repacking (default=-1).' )
parser.add_option('--neighbor', dest = 'neighbor',
    default = -1,    # default to False
    help = '>0 --> CB-CB distance to define loop neighborhood for repacking only, 0 --> loop-only, <0 --> repack neighborhood equal to minimization one (default=-1).' )
parser.add_option('--minneighbor', dest = 'minneighbor',
    default = 0,    # default to False
    help = 'CB-CB distance (loop-environment) to define loop neighborhood for MinMover descent minimization, set 0 to loop-only (default=0).' )
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
    default = 'ref2015',    # default energy function
    help = 'Rosetta energy function, e.g: score12, talaris2013, talaris2014, or ref2015' )
parser.add_option('--verbose', dest = 'verb',
    default = 'error',    # default verbose level
    help = 'Verbose levels are fatal/error/warning/info/debug/trace, (default = error)' )

parser.add_option('--refine', dest = 'refine',
    default = '',    # default to False
    help = 'Enables a Rosetta\'s KIC (very slow) or CCD refinements (slow) (default = disabled)' )
parser.add_option('--nbest', dest = 'nbest',
    default = '',    # default to False
    help = 'Number of best solutions to store' )
parser.add_option('--nrep', dest = 'nrep',
    default = 1,    # default to False
    help = 'Number of minimization/refinement repetitions to find the lowest energy and store the corresponding loop' )
parser.add_option('--refine_best', dest = 'refine_best',
    default = '',    # default to False
    help = 'Enables: fast gradient-minimization (monrelax), old gradient-minimization (fastrelax), very slow Rosetta\'s KIC (KIC), or slow CCD (CCD) refinements for best solutions (default = disabled)' )
parser.add_option('--nmin_best', dest = 'nmin_best',
    default = 5,    # default to False
    help = 'Number of minimization steps of the Monte-Carlo plus Minimization for best loops (only for \"fastrelax\") (default=5)' )
parser.add_option('--fas', dest = 'fas',
    default = '',    # default to False
    help = 'VdW repulsive factors for fast-relax minimization protocol (default= 0.05,0.30,0.60,1.00)' )
parser.add_option('--tols', dest = 'tols',
    default = '',    # default to False
    help = 'Convergence/tolerance factors for quick-relax minimization protocol (default= 0.01,0.01,0.001,0.00001)' )
parser.add_option('--server', dest="server", 
    default=False, 
    action="store_true",
    help = 'Enable RCD-server output:\n\t <base>_dhs_ros_##.txt (refined dihedral angles)\t\t\t(default=disabled)' )
parser.add_option('--ntpsi', dest="ntpsi", 
    default=False, 
    action="store_true",
    help = 'Enable N-terminal-PSI dihedral angle optimization. (default=disabled)' )
parser.add_option('--native', dest="native", 
    default=False, 
    action="store_true",
    help = 'Computes the neighborhood mask from the native loop (that must be present in the input PDB). Otherwise, the mask will be computed as those residues whose CB atom is closer than --minneighbor distance from any loop CB. (default=disabled)' )
parser.add_option('--nminneigh', dest = 'nminneigh',
    default = '',    # default to False
    help = 'Number of loops used to define --minneighbor neighborhood. By default, all loops will be screened. (default=all)' )

parser.add_option('--premature', dest="premature", 
    default=False, 
    action="store_true",
    help = 'Forces premature exit from outher for-loops in gradient-minimizer \"monrelax\". (default=disabled)' )
parser.add_option('--precursor', dest="precursor", 
    default=False, 
    action="store_true",
    help = 'Uses the \"precursor\" structure instead of the refined one. ONLY for \"best\" protocol. (default=disabled)' )


(options,args) = parser.parse_args()

# os.environ['LD_LIBRARY_PATH']='/home/mon/prgs/PyRosetta.Ubuntu-12.04LTS.64Bit.monolith/'

import time
start = time.time()
import os;

# PyRosetta-4 
import pyrosetta
import pyrosetta.rosetta as rosetta
from pyrosetta import init, PyMOLMover, Pose, pose_from_file, ScoreFunction, create_score_function, get_fa_scorefxn, Vector1
from pyrosetta.rosetta import core, protocols

# Import my stuff into current namespace...
from monfuncs import *

# WARNING: option '-constant_seed' is for testing only! MAKE SURE TO REMOVE IT IN PRODUCTION
#init(extra_options = "-constant_seed")

verb = options.verb
if not verb:
	verb = 'error'
print("Using %s verbose level." % verb)

energy = options.energy
if not energy:
	energy = 'ref2015'
print("Using Rosetta's %s energy function (for repacking and minimization)." % energy)
talaris = ""
if "talaris" in energy: 
	talaris = " -corrections::restore_talaris_behavior" # required to work with Talaris energies
	print("...adding \"%s\" option to PyRosetta init." % talaris)

init(extra_options = "-ex2 -rebuild_disulf false -detect_disulf false -ignore_zero_occupancy false -allow_omega_move false -ignore_unrecognized_res -out:levels core:" + verb + " -out:levels protocols:" + verb + talaris)
# init(extra_options = "-ex2 -rebuild_disulf false -detect_disulf false -ignore_zero_occupancy false -allow_omega_move false -ignore_unrecognized_res -corrections::restore_talaris_behavior -out:levels core:" + verb + " -out:levels protocols:" + verb)
#
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
#
# About tracers verbose: https://www.rosettacommons.org/node/10018
# the "out" (verbose) levels are fatal/error/warning/info/debug/trace with the default being 'info'


end = time.time()
print("Init time: %f" % (end - start))

# Loading parser variables

base_in = options.base_in

pdb_in = options.pdb_in
if not pdb_in:
	pdb_in = "%s.pdb" % base_in # Take input pdb file name from input base-name
if size > 1:
	base_out = "%s_%d" % (options.base_out,rank+1)  # Each worker will produce its own output...
else:
	base_out = options.base_out
print("OUTPUT BASENAME: %s\n" % base_out)

if rank == 0:	# Only master process should write non-output files
	pdb_read = "%s_input.pdb" % base_out # readed pdb with all atoms added by PyRosetta

dh_file = options.dh_file
if not dh_file:
	dh_file = "%s_dh.txt" % base_in # Take input dihedrals file name from input base-name

val_file = options.val_file
if not val_file:
	val_file = "%s_val.txt" % base_in # Take input valence angles file name from input base-name

len_file = options.len_file
if not len_file:
	len_file = "%s_len.txt" % base_in # Take input bond lengths file name from input base-name

ri = int(options.ri)
rf = int(options.rf)
chain = options.chain

rmsds_file = options.rmsds_file
if not rmsds_file:
	rmsds_file = "%s_rmsd.txt" % base_in # Take input rmsds file name from input base-name

Bfactor = options.Bfactor

nmin = int(options.nmin)

nmin_best = int(options.nmin_best)

loopI = options.loopI
loopF = options.loopF
if not loopI:
	loopI = 1
else:
	loopI = int(loopI)

normsds = options.normsds
if not normsds:
	normsds = 0
else:
	normsds = 1


initneighbor = options.initneighbor
if not initneighbor:
	initneighbor = int(-1)
else:
	initneighbor = float(initneighbor)
print("Parser: initneighbor= %f" % initneighbor)

neighbor = options.neighbor
if not neighbor:
	neighbor = int(-1)
else:
	neighbor = float(neighbor)
print("Parser: neighbor= %f" % neighbor)

minneighbor = options.minneighbor
if not minneighbor:
	minneighbor = int(0)
else:
	minneighbor = float(minneighbor)
print("Parser: minneighbor= %f" % minneighbor)

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
	print("Parser: Loop-Repacking disabled (irepack=%d)" % irepack)
	irepack = int(-1); # Repacking disabled

refine = options.refine

nbest = options.nbest
if not nbest:
	nbest = int(0)
else:
	nbest = int(options.nbest)
print("Parser: nbest= %d" % nbest)

refine_best = options.refine_best

nrep = options.nrep
if not nrep:
	nrep = int(1)
else:
	nrep = int(options.nrep)

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
print("Using %s flavor in descent minimization (mintype option for MinMover)." % mintype)

server = options.server
print('SERVER OUTPUT: ', server)
ntpsi = options.ntpsi
print('N-terminal PSI angle optimization: ', ntpsi)

native = options.native
if not native:
	print('PARSER: Neighborhood mask computed as those residues whose CB atom is closer than --minneighbor distance from any loop CB.')
else:
	print('PARSER: Neighborhood mask taken from the native loop (that must be already present in the input PDB).')

energy_file = "%s_rosetta.txt" % base_out # Output table data (from Rosetta)

premature = options.premature
print('Parser: premature= ', premature)
precursor = options.precursor
print('Parser: precursor= ', precursor)


#######################################
# MAIN
#######################################

# Load the data from pdb_file into the pose
print("Loading input pdb: %s" % pdb_in)
start = time.time()
pose = Pose() # Create an empty Pose object
pose = pose_from_file(pdb_in) # Warning, this automatically rebuilds missing atoms and repacks...
end = time.time()
print("Load time: %f" % (end - start))

# Show pdb info
info = pose.pdb_info()
print(info)

# write input pdb (now full atom and packed)
if size > 1:
	if server and rank == 0: # Only for master MPI process
		pose.dump_pdb(pdb_read)
else:
	if server and loopI == 1: # Only for the first thread (if crap-multi-threading enabled in server)
		pose.dump_pdb(pdb_read)

# Find chain ID if not provided by parser
chain_num = 1
if not chain or len(chain) != 1:
	chain = info.chain(chain_num) # Get chain ID of first residue
	print("Chain %c selected by default" % chain)
else:
	print("Chain %c selected" % chain)

# Find correspondence between PDB and Rosetta numbering schemes for residues
rri = info.pdb2pose(chain,ri) # initial Rosetta residue index
if rri == 0:
	print("Residue %d from chain %c not found! Forcing exit!" % (ri, chain))
	exit()
rrf = info.pdb2pose(chain,rf) # final Rosetta residue index
if rrf == 0:
	print("Residue %d from chain %c not found! Forcing exit!" % (rf, chain))
	exit()

if (rrf-rri) != (rf-ri):
	print("ERROR: Correspondence mismatch between PDB and Rosetta numbering! (rrf=%d rri=%d rf=%d ri=%d)" % (rrf,rri,rf,ri))
	exit()

print("PDB numbered segment %d-%d corresponds to %d-%d in Rosetta." % (ri, rf, rri, rrf))


loop_sel = range( rri+1 , rrf ) # loop range selection (for RMSD computations) NOT-INCLUDING ANCHORS
print("Selected loop range from %d to %d (residues considered in RMSDs)" % (rri+1,rrf-1))
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
BBCB_rmsd = []              ## Added 2023 (ivan) RMSD like 'BBOCB_rmsd' whitout Oxigen (implemented in all outputs below)
BBCB_rmsd.append('N')
BBCB_rmsd.append('CA')
BBCB_rmsd.append('C')
BBCB_rmsd.append('CB')
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


# NATRO # use native amino acid and native rotamer (does not repack)
# NATAA # use native amino acid but allow repacking to other rotamers
# PIKAA ILV # use only the following amino acids and allow repacking between them
# ALLAA # use all amino acids and all repacking
print("Generating resfile from PDB numbered segment")
specific_design = {}
for r in range( rri+1 , rrf ): # Rosetta numbering
     specific_design[r] = 'NATAA' # NATAA use native amino acid but allow repacking to other rotamers
#     specific_design[r] = 'ALLAA' # ALLAA use all amino acids and all repacking

# create a standard ScoreFunction
# scorefxn = get_fa_scorefxn('talaris2013') #  create_score_function_ws_patch('standard', 'score12')
# scorefxn = create_score_function_ws_patch('standard', 'score12')

# scorefxn = get_fa_scorefxn() #  create_score_function_ws_patch('standard', 'score12')
# scorefxn = create_score_function('talaris2013_calibrated')
#scorefxn = create_score_function('talaris2013') # state-of-the-art Rosetta's energy
#scorefxn = create_score_function('talaris2014') # state-of-the-art Rosetta's energy
# scorefxn = create_score_function('cacagorda') # state-of-the-art Rosetta's energy (this does not work...)
scorefxn = create_score_function(energy) # state-of-the-art Rosetta's energy

#scorefxn.apply_patch_from_file("score12")

# MON FIXED (12/1/2016)
scorefxn.set_weight(core.scoring.chainbreak,100.0) # 100 (instead of 1) is important to avoid chain breaks
scorefxn.set_weight(core.scoring.coordinate_constraint,1.0)

print(scorefxn)
# print "the fa_atr= %f" % scorefxn[fa_atr]

# write a resfile to select loop aminoacids for repacking
#resfile_name = "%s_resfile.txt" % (base_out)
#generate_resfile_from_pose(pose, resfile_name, False, specific = specific_design)

# Compute native loop Energy (initial pose should have native loop)
poseE0 = scorefxn(pose)
poseE_res = [] # array with residue energies
for i in range(rri,rrf+1):
	poseE_res.append(pose.energies().residue_total_energies(i).sum() / Bfactor)

# This trick backups the loop PDBinfo from input molecule
dummy = Pose()
dummy.assign(pose) # full copy of "pose"

# Deleting non-loop (and non-anchor) residues 
# (this seems the only way to assing later a valid PDBinfo object to each loop...)
print("Deleting residues from %d to %d in dummy" % (rrf+1,pose.total_residue()))
rosetta.protocols.grafting.delete_region(dummy,rrf+1,pose.total_residue())
#dummy.conformation().delete_residue_range_slow(rrf+1,pose.total_residue())
#dummy.pdb_info().obsolete(false)
print("Deleting residues from %d to %d in dummy" % (1,rri-1))
rosetta.protocols.grafting.delete_region(dummy,1,rri-1)
#dummy.conformation().delete_residue_range_slow(1,rri-1)
#dummy.pdb_info().obsolete(false)

#print "Retrieve dummy loop: residues from %d to %d" % (rri,rrf)
#dummy = rosetta.protocols.grafting.return_region(pose, rri, rrf)

#dummy.conformation().detect_disulfides() # This is mandatory, since delete_region seems not to refresh disulfides...

#print "Chain %c" % info.chain(dummy.residue(1).chain())

loop_info = dummy.pdb_info() # Store valid PDBinfo to be asigned later to each generated loop

# Apply PDBinfo from dummy...
#dummy.pdb_info(loop_info)

# Save native loop
if server and loopI == 1: # Only for the first thread (if crap-multi-threading enabled in server)
	if rank == 0:	# Only master process should write non-output files
		loop_name = "%s_native.pdb" % (base_out)
		write_pdb(dummy, loop_name, poseE_res)


# Native Pose backup (to compute RMSDs)
posen = Pose()
posen.assign(pose) # Create a full copy the initial pose
posen.pdb_info(info)

loop_name = "%s_native.pdb" % (base_out)
#write_pdb(dummy, loop_name, poseE_res)

# INITIAL STUFF...
if initneighbor >= 0:
	print("Repacking using Initial neighborhood (initneighbor=%f)\n" % (initneighbor))
	initneighbors = neighbors(pose,initneighbor,rri,rrf)  # Define initial neighborhood for initial pose
	print(initneighbors)
	neighbor_packer(pose,scorefxn,initneighbors)  # Repack current pose
	print("Minimization with Initial neighborhood (initneighbor=%f)\n" % (initneighbor))
	monrelax(pose,scorefxn,nmin,-1,initneighbors,fas,tols,rri,rrf,repackouther,irepack,ntpsi,mintype,premature) # "nmin" iters with all "fas" (slower)
	# Write non-repacked pdb
	if rank == 0:	# Only master process should write non-output files
		initmin_name = "%s_initmin.pdb" % (base_out)
		pose.dump_pdb(initmin_name)


# Compute Energy of the Minimized initial structure (initial pose should have native loop or not)
poseEmin = scorefxn(pose)

# Compute RMSDs of the Minimized initial structure
if normsds == 0:
	ca_rmsd = rmsd(posen,pose,loop_sel,CA_rmsd)
	bb_rmsd = rmsd(posen,pose,loop_sel,BB_rmsd)
	bbo_rmsd = rmsd(posen,pose,loop_sel,BBO_rmsd)
	bbcb_rmsd = rmsd(posen,pose,loop_sel,BBCB_rmsd)
	bbocb_rmsd = rmsd(posen,pose,loop_sel,BBOCB_rmsd)
	ha_rmsd = rmsd(posen,pose,loop_sel,None,HA_rmsd)
	sc_rmsd = rmsd(posen,pose,loop_sel,None,SC_rmsd)
	full_rmsd = rmsd(posen,pose,loop_sel,None,None)

# Reading Internal Coordinates
print("Reading Internal Coordinates from %s, %s, and %s files." % (dh_file,val_file,len_file))
dihedrals = read_dihedrals(dh_file) # Read dihedral angles from file
valences = read_valences(val_file) # Read valence angles from file
lengths = read_lengths(len_file) # Read bond lengths from file

# This must be done after reading "dihedrals" file
if not loopF:
	loopF = int(len(dihedrals)) # Loops
else:
	loopF = int(loopF)

# Load per thread is computed here, it must be done after reading "dihedrals" file
if size > 1: # Parallel stuff...
	lpc = int(loopF/size) # Loops per core
	r = loopF%size # Remainder
	if rank+1 <= r:
		loopI = rank * lpc + rank + 1
		loopF = loopI + lpc
	else:
		loopI = rank * lpc + r + 1
		loopF = loopI + lpc - 1

print("Process rank %3d loop range: %d to %d\n" % (rank,loopI,loopF))


# Reading RMSDs file
if rmsds_file:
	print("Reading RMSDs and KORP energies file %s" % (rmsds_file))
	rmsds = read_column(rmsds_file,1,'#')
	korps = read_column(rmsds_file,5,'#')


if nbest <= 0:
	# Deleting previous loops (Multi-PDB file)
	loops_name = "%s_loops.pdb" % (base_out)
	print("Deleting previous loops (Multi-PDB file): %s" % (loops_name))
	f = open(loops_name,'w') # Deleting previous loops (Multi-PDB file)
	f.close()
	# Chapa to get KORP's native energy...
	b = open(rmsds_file, "r").readlines()
	line = b[1].split() #Gets the second last value from the list and split on whitespace
	korpEnat = float(line[6])
	# b.close()
	print("KORP's native energy: %f" % (korpEnat))
	# Open output file with energies
	f = open(energy_file, 'w+')
	f.write("#%4s %8s %8s %6s %12s %5s %5s %5s %8s %5s %5s %5s %5s %5s %5s %5s %5s\n" % ('loop', 'E_full', 'E_loop', 'R_rcd', 'E_KORP', 'R_BB', 'R_HA', 'R_SC', 'E_full2', 'R_CA', 'R_BB', 'R_BBO','R_CB', 'R_OCB', 'R_HA', 'R_SC', 'R_ALL') )
	f.write("#%4s %8.2f %8.2f %6.3f %12f %5.2f %5.2f %5.2f %8.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n" % (' ', poseE0, scorefxn(dummy), 0.0, korpEnat, 0.0, 0.0, 0.0, poseEmin, ca_rmsd, bb_rmsd, bbo_rmsd, bbcb_rmsd, bbocb_rmsd, ha_rmsd, sc_rmsd, full_rmsd) )

# Storing dihedral angles of the initial loop (to get anchors dihedrals) in memory (Server output)
if server:
	dihedrals2 = [] # Array to store refined dihedral angles
	dihedrals0 = get_dihedrals(pose,rri,rf-ri+1)

best_list = [] # Initialize Best poses list

# Minimization neighborhood definition...
if native:
	# Define CONSTANT neighborhood for ALL loops minimization
	print("NEIGHBORHOOD CALCULATION FROM INITIAL LOOP (input PDB loop)")
	minneighbors = neighbors(pose,minneighbor,rri,rrf)
else: 
	nminneigh = options.nminneigh
	if not nminneigh:
		nminneigh = int(len(dihedrals))
	else:
		nminneigh = int(nminneigh)
		if nminneigh > int(len(dihedrals)):
			nminneigh = int(len(dihedrals))		

	# Selection (the best choice)
	print("NEIGHBORHOOD CALCULATION FROM the first %d SAMPLED LOOPS" % nminneigh)
	start = time.time()
	# Create a full copy the initial pose for every loop (Avoid repacking errors between different loops)
	pose2 = Pose()
	pose2.assign(pose)
	minneighbors = 0  # Neighbors' hash is first initializaed with 0
	for l in range(0,nminneigh):

		# Update Internal Coordinates from loaded ICs arrays
		update_dihedrals(pose2, dihedrals[l], rri)
		update_valences(pose2, valences[l], rri)
		update_lengths(pose2, lengths[l], rri)

		# Neighbors' hash is incrementally built...
		minneighbors = neighbors(pose2,minneighbor,rri,rrf,minneighbors)
		print("Neighborhood from loop %d = %d resiudes\n" % (l+1,len(minneighbors)))
	end = time.time()
	print("Initial Neighborhood detection time: %f s" % (end-start))
	

#myneighbors = neighbors(pose,10,rri,rrf)  # Define neighborhood for current pose

# UPDATE INTERNAL COORDINATES
#  pose --> PyRosetta's pose
#  pdbs --> Cartesian coordinates 4D array: pdbs[loop][atom][residue][coord]
for l in range(loopI-1,loopF):
	start = time.time()
	print("Processing loop %d\n" % l)

	# Create a full copy the initial pose for every loop (Avoid repacking errors between different loops)
	pose2 = Pose()
	pose2.assign(pose)

	# Update Internal Coordinates from loaded ICs arrays
	update_dihedrals(pose2, dihedrals[l], rri)
	update_valences(pose2, valences[l], rri)
	update_lengths(pose2, lengths[l], rri)

	# print "Dihedrals as loaded from IC files..."
	# for q in range(0,rrf-rri+1): # iters from resi to resi+nres residues
	#    print "%12f   %12f" % (pose2.phi(rri+q),pose2.psi(rri+q))  # Phi/Psi angles

	# Write non-repacked pdb
	# out_name = "pre_%05d.pdb" % (l+1)
	# pose2.dump_pdb(out_name)

	#loop_name = "repack_%05d.pdb" % (l+1)
	#pose2.dump_pdb(loop_name)

	# Computing Rosetta Energy (full protein)
	#fullE = scorefxn(pose2)
	bb_rmsd0 = 0.0
	ha_rmsd0 = 0.0
	sc_rmsd0 = 0.0
	ca_rmsd = 0.0
	bb_rmsd = 0.0
	bbo_rmsd = 0.0
	bbcb_rmsd = 0.0
	bbocb_rmsd = 0.0
	ha_rmsd = 0.0
	sc_rmsd = 0.0
	full_rmsd = 0.0

	if normsds == 0:
		bb_rmsd0 = rmsd(posen,pose2,loop_sel,BB_rmsd)
		ha_rmsd0 = rmsd(posen,pose2,loop_sel,None,HA_rmsd)
		sc_rmsd0 = rmsd(posen,pose2,loop_sel,None,SC_rmsd)
	mypose = pose2 # Some trick...

#	if nmin <= 0: # If no minimization "nmin<=0" (only repacking)
#		fullE2 = scorefxn(pose2)
#	else:
#		fullE2 = 1e10 # some very-high number

	fullE2 = scorefxn(pose2)
	fullE = fullE2

	# Repeat minimization/refinement protocol "nrep" times retaining the lowest energy loop
	for r in range(1,nrep+1):
		start2 = time.time()

		# Hard copy of initial pose
		posel = Pose()
		posel.assign(mypose) 

		# MON TEST (12/1/2016)
		# Loop must be repacked always (otherwise cheating...)
#		print "Initial Repacking of loop %d (repeat %d) only-loop (neighbor=0)\n" % (l,r)
#		myneighbors = neighbors(posel,50,rri,rrf)  # Define neighborhood for current pose
#		neighbor_packer(posel,scorefxn,myneighbors)  # Repack current pose
		
		# Quick-relax stepest descent protocol
		if nmin > 0:
			# print "rri= %d  rrf= %d" % (rri,rrf)
			monrelax(posel,scorefxn,nmin,neighbor,minneighbors,fas,tols,rri,rrf,repackouther,irepack,ntpsi,mintype,premature)

			# Print Phi/Psi angles
			# for q in range(0,rrf-rri+1): # iters from resi to resi+nres residues
			#	print "%12f %12f   %12f %12f" % (pose2.phi(rri+q),posel.phi(rri+q),pose2.psi(rri+q),posel.psi(rri+q))  

		# END of quick-relax stepest descent protocol
		
		# KIC / CCD final refinement protocol
		if refine:
			loop = Loop(rri,rrf)
			loop.auto_choose_cutpoint(posel)
			set_single_loop_fold_tree(posel,loop)
			add_cutpoint_variants(posel)
			loops = Loops()
			loops.add_loop(loop)
			loop_refine = 0
			if 'KIC' in refine:
				loop_refine = rosetta.protocols.loops.loop_mover.refine.LoopMover_Refine_KIC( loops )
			if 'CCD' in refine:
				loop_refine = rosetta.protocols.loops.loop_mover.refine.LoopMover_Refine_CCD( loops )
				loop_refine.temp_initial(0.6)
				loop_refine.temp_final(0.1)
				loop_refine.outer_cycles(5)
				loop_refine.max_inner_cycles(5)
			loop_refine.apply( posel )
		# END of KIC / CCD final refinement protocol

		myenergy = scorefxn(posel)
		if fullE2 > myenergy:
			pose2 = posel  # pose2 updaterd by the best minimized/refined pose
			fullE2 = myenergy  # update best minimized/refined energy
		end2 = time.time()
		print("*** Repetition %2d energy: %f,  best energy: %f (%f s)" % (r, myenergy, fullE2, end2-start2))

	resE = [] # array with residue energies
	for i in range(rri,rrf+1):
		resE.append(pose2.energies().residue_total_energies(i).sum() / Bfactor)
	# Apply PDBinfo from pdb...
	
	
	if normsds == 0:
		ca_rmsd = rmsd(posen,pose2,loop_sel,CA_rmsd)
		bb_rmsd = rmsd(posen,pose2,loop_sel,BB_rmsd)
		bbo_rmsd = rmsd(posen,pose2,loop_sel,BBO_rmsd)
		bbcb_rmsd = rmsd(posen,pose2,loop_sel,BBCB_rmsd)
		bbocb_rmsd = rmsd(posen,pose2,loop_sel,BBOCB_rmsd)
		ha_rmsd = rmsd(posen,pose2,loop_sel,None,HA_rmsd)
		sc_rmsd = rmsd(posen,pose2,loop_sel,None,SC_rmsd)
		full_rmsd = rmsd(posen,pose2,loop_sel,None,None)

	# Trim loop region from full-protein pose
#	loop = rosetta.protocols.grafting.return_region(pose2,rri,rrf)
#	loop.conformation().detect_disulfides() # This is mandatory, since delete_region seems not to refresh disulfides...

	# Mon: Modified 13/10/2015
	# Deleting non-loop residues 
	# (this seems the only way to assign later a valid PDBinfo object to each loop...)
	loop = Pose()
	loop.assign(pose2);
	print("Deleting residues from %d to %d in loop" % (rrf+1,loop.total_residue()))
	rosetta.protocols.grafting.delete_region(loop,rrf+1,loop.total_residue())
	print("Deleting residues from %d to %d in loop" % (1,rri-1))
	rosetta.protocols.grafting.delete_region(loop,1,rri-1)

	# Apply PDBinfo from dummy...
	loop.pdb_info(loop_info)
	



	# Compute loop Energy
	#loopE = 0.0
	loopE = scorefxn(loop) # loop energy should be updated

	# Computing per-residue energy
	end = time.time()
	print("Relaxed Loop %5d energy: E= %f  E2= %f  (%f s)" % (l+1, fullE, fullE2, end-start))

	# Write repacked pdb (non-trimmed)
	#out_name = "repack_%04d.pdb" % (l+1)
	#pose2.dump_pdb(out_name)

	# Save PDB
	#loop_name = "%s_%05d.pdb" % (base_out, l+1)
	#print "KK: %s\n" % (loop_name)
	#loop.dump_pdb(loop_name)
	#write_pdb(loop, loop_name, resE)
	if nbest <= 0:
		print("Storing loop %d, residues from %d to %d" % (l+1, rri, rrf))
		append_mpdb(loop, loops_name, resE, l+1)

	# Storing dihedral angles of refined loops in memory (Server output)
	if server:
		dihedrals2.append( get_dihedrals(loop) )

	myrmsd = 0.0
	mykorp = 0.0
	if rmsds_file:
		myrmsd = float(rmsds[l])
		mykorp = float(korps[l])

	if precursor and nbest > 0:
		pose2 = mypose # pose2 is back...

	# Keep RMSDs....	
	data = "%5d %8.2f %8.2f %6.3f %12f %5.2f %5.2f %5.2f %8.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f" % (l+1, fullE, loopE, myrmsd, mykorp, bb_rmsd0, ha_rmsd0, sc_rmsd0, fullE2, ca_rmsd, bb_rmsd, bbo_rmsd, bbcb_rmsd, bbocb_rmsd, ha_rmsd, sc_rmsd, full_rmsd)

	if nbest > 0:
#		data = "%5d %8.2f %8.2f %6.3f %5.2f %5.2f %5.2f %8.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f" % (l+1, fullE, loopE, myrmsd, bb_rmsd0, ha_rmsd0, sc_rmsd0, fullE2, ca_rmsd, bb_rmsd, bbo_rmsd, bbocb_rmsd, ha_rmsd, sc_rmsd, full_rmsd)
		# satan = scorefxn(pose2)  # Initial refined energy
		# print "loop %d --> energy= %f" % (l+1,satan)
		keepbest(nbest,best_list,pose2,fullE2,data)
	else:
#		f.write("%5d %8.2f %8.2f %6.3f %5.2f %5.2f %5.2f %8.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f\n" % (l+1, fullE, loopE, myrmsd, bb_rmsd0, ha_rmsd0, sc_rmsd0, fullE2, ca_rmsd, bb_rmsd, bbo_rmsd, bbocb_rmsd, ha_rmsd, sc_rmsd, full_rmsd))
		f.write(data+"\n")
		f.flush()


# Saving dihedral angles of refined loops as plain text format (Server output)
if server:
	for j in range(0,len(dihedrals2[0])//2): # Screen residues # "//" floor division operator in Python>3.x
		refdhs_name = "%s_refdhs_%02d.txt" % (base_out,j+1)
		f=open(refdhs_name,'w') # Deleting previous loops (Multi-PDB file)
		if j == 0: # Nt-anchor
			for i in range(0,len(dihedrals2)): # Screen refined loops
				f.write("%10.5f %10.5f\n" % (dihedrals0[j*2],dihedrals2[i][j*2+1]) )
		elif j == rf-ri: # Ct-anchor
			for i in range(0,len(dihedrals2)): # Screen refined loops
				f.write("%10.5f %10.5f\n" % (dihedrals2[i][j*2],dihedrals0[j*2+1]) )
		else: # Non-anchors
			for i in range(0,len(dihedrals2)): # Screen refined loops
				f.write("%10.5f %10.5f\n" % (dihedrals2[i][j*2],dihedrals2[i][j*2+1]) )
		f.close()

# "BEST" protocol...
if nbest > 0:
	# loops_name = "%s_best.pdb" % (base_out)
	loops_name = "%s_loops.pdb" % (base_out)
	print("Deleting previous loops (Multi-PDB file): %s" % (loops_name))
	f=open(loops_name,'w') # Deleting previous loops (Multi-PDB file)
	f.close()

	# Sorting best solutions
	best_list.sort(key=getKey, reverse=False)  # Sort list
	# print best_list

	# Final refinement with "best" candidates
	for i in range(0,len(best_list)):
		start = time.time()
		if refine_best:

			if 'monrelax' in refine_best:
				# monrelax(best_list[i][1],scorefxn,nmin,neighbor,minneighbors,fas,tols,rri,rrf,repackouther,irepack,ntpsi,mintype,premature)
				monrelax(best_list[i][1],scorefxn,nmin_best,neighbor,minneighbors,fas,tols,rri,rrf,False,0,ntpsi,mintype,False)
			else:
				# mierda0 = scorefxn(best_list[i][1])  # Initial refined energy

				loop = Loop(rri,rrf)
				loop.auto_choose_cutpoint(best_list[i][1])
				set_single_loop_fold_tree(best_list[i][1],loop)
				add_cutpoint_variants(best_list[i][1])
				loops = Loops()
				loops.add_loop(loop)
				loop_refine = 0

				# mierda = scorefxn(best_list[i][1])  # Initial refined energy

				if 'fastrelax' in refine_best:
					pd2relax3(best_list[i][1],scorefxn,nmin_best,neighbor)
				if 'KIC' in refine_best:
					loop_refine = rosetta.protocols.loops.loop_mover.refine.LoopMover_Refine_KIC( loops )
					loop_refine.apply( best_list[i][1] ) 
				if 'CCD' in refine_best:
					loop_refine = rosetta.protocols.loops.loop_mover.refine.LoopMover_Refine_CCD( loops )
					loop_refine.temp_initial(0.6)
					loop_refine.temp_final(0.1)
					loop_refine.outer_cycles(5)
					loop_refine.max_inner_cycles(5)
					loop_refine.apply( best_list[i][1] ) 

			best_list[i][3] = scorefxn(best_list[i][1])  # Final refined energy

		if normsds == 0:
			ca_rmsdF = rmsd(pose,best_list[i][1],loop_sel,CA_rmsd)
			bb_rmsdF = rmsd(pose,best_list[i][1],loop_sel,BB_rmsd)
			bbo_rmsdF = rmsd(pose,best_list[i][1],loop_sel,BBO_rmsd)
			bbcb_rmsdF = rmsd(pose,best_list[i][1],loop_sel,BBCB_rmsd)
			bbocb_rmsdF = rmsd(pose,best_list[i][1],loop_sel,BBOCB_rmsd)
			ha_rmsdF = rmsd(pose,best_list[i][1],loop_sel,None,HA_rmsd)
			sc_rmsdF = rmsd(pose,best_list[i][1],loop_sel,None,SC_rmsd)
			full_rmsdF = rmsd(pose,best_list[i][1],loop_sel,None,None)
			data = "%8.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f" % (best_list[i][3], ca_rmsdF, bb_rmsdF, bbo_rmsdF, bbcb_rmsdF, bbocb_rmsdF, ha_rmsdF, sc_rmsdF, full_rmsdF)
			best_list[i][4] = data
		
		end = time.time()
#		print "Best loop %d refined in %.3f s, from energy %.3f (%.3f %.3f) to %.3f, RMSDs %s\n" % (i,end-start,best_list[i][0],mierda0,mierda,best_list[i][3],data)
		print("Best loop %d refined in %.3f s, from energy %.3f to %.3f, RMSDs %s\n") % (i,end-start,best_list[i][0],best_list[i][3],data)

	# Sorting best solutions again
	best_list.sort(key=getKey2, reverse=False)  # Sort list

	for i in range(0,len(best_list)):
		# Get individual residue energies
		resE = [] # array with residue energies
		for j in range(rri,rrf+1):
			resE.append(best_list[i][1].energies().residue_total_energies(j).sum() / Bfactor)

		# Trim loop region from full-protein pose
#		loop = rosetta.protocols.grafting.return_region(best_list[i][1],rri,rrf)
#		loop.conformation().detect_disulfides() # This is mandatory, since delete_region seems not to refresh disulfides...

		# Mon: Modified 13/10/2015
		# Deleting non-loop residues 
		# (this seems the only way to assing later a valid PDBinfo object to each loop...)
		loop = Pose()
		loop.assign(best_list[i][1]);
		print("Deleting residues from %d to %d in loop") % (rrf+1,loop.total_residue())
		rosetta.protocols.grafting.delete_region(loop,rrf+1,loop.total_residue())
		print("Deleting residues from %d to %d in loop") % (1,rri-1)
		rosetta.protocols.grafting.delete_region(loop,1,rri-1)

		# Apply PDBinfo from dummy...
		loop.pdb_info(loop_info)
			
		# Get refined loop index to be used as Multi-PDB index. (dirty but works...)
		myindex = int(best_list[i][2].split()[0])
		print("myindex= %d") % myindex

		# Save Multi-PDB
		append_mpdb(loop, loops_name, resE, myindex)
		print("Storing best loop %d, residues from %d to %d") % (i+1, rri, rrf)

	# best_name = "%s_best.txt" % (base_out)
	best_name = energy_file
	print("Saving best loops data file: %s") % (best_name)
	b=open(best_name,'w') 
	b.write("#%4s %8s %8s %6s %5s %5s %5s %8s %5s %5s %5s %5s %5s %5s %5s %5s " % ('loop', 'E_full', 'E_loop', 'R_rcd', 'R_BB', 'R_HA', 'R_SC', 'E_full2', 'R_CA', 'R_BB', 'R_BBO', 'R_CB', 'R_OCB', 'R_HA', 'R_SC', 'R_ALL') )
	b.write("%8s %5s %5s %5s %5s %5s %5s %5s %5s\n" % ('E_full3', 'R_CA', 'R_BB', 'R_BBO', 'R_CB', 'R_OCB', 'R_HA', 'R_SC', 'R_ALL') )
	for i in range(0,len(best_list)):
		b.write("%s %s\n" % (best_list[i][2],best_list[i][4]))
	b.close()


exit()


#########################################################
# ENERGY STUFF...
#########################################################

start = time.time()

# Get H-bonding data...
#hbond_set = rosetta.core.scoring.hbonds.HBondSet()
#pose.update_residue_neighbors()
#rosetta.core.scoring.hbonds.fill_hbond_set(pose,False,hbond_set)
#hbond_set.show(pose)
#print hbond_set
#from toolbox import get_hbonds
#hset = get_hbonds(pose)
# hset.show(pose) # Shows all H-bonds
# hset.show(pose, 97)
# res97 = hset(pose, 97)
# print "res97 energy= %f" % res97[energy]

# USE LATER: To store text output in string
# b = rosetta.utility.OStringStream()
# sf.show(b, pose)
# print b.str()

#scorefxn = get_fa_scorefxn()
print("%4s %11s %11s %11s %12s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s %11s") % ('#res','fa_atr', 'fa_rep', 'fa_sol', 'fa_intra_rep', 'fa_elec', 'fa_dun', 'dslf_fa13', 'rama', 'omega', 'p_aa_pp', 'hbond_sr_bb', 'hbond_lr_bb', 'hbond_bb_sc', 'hbond_sc', 'pro_close', 'ref')

res_energy = 0.0 # Total energy summed over all individual terms from all residues

tot_fa_atr = 0.0 # 
tot_fa_rep = 0.0 # 
tot_fa_sol = 0.0 # 
tot_fa_intra_rep = 0.0 # 
tot_fa_elec = 0.0 # 
tot_fa_dun = 0.0 # 
tot_dslf_fa13 = 0.0 # 
tot_rama = 0.0 # 
tot_omega = 0.0 # 
tot_p_aa_pp = 0.0 # 
tot_hbond_sr_bb = 0.0 # 
tot_hbond_lr_bb = 0.0 # 
tot_hbond_bb_sc = 0.0 # 
tot_hbond_sc = 0.0 # 
tot_pro_close = 0.0 # 
tot_ref = 0.0 # 


for i in range(1,pose.total_residue()+1):
# for i in range( rri , rrf+1 ):
#	pose.energies().residue_total_energies(res_num)[fa_rep]
	resE = pose.energies().residue_total_energies(i)
	tot_fa_atr += resE[fa_atr]
	tot_fa_rep += resE[fa_rep]
	tot_fa_sol += resE[fa_sol]
	tot_fa_intra_rep += resE[fa_intra_rep]
	tot_fa_elec += resE[fa_elec]
	tot_fa_dun += resE[fa_dun]
	tot_dslf_fa13 += resE[dslf_fa13]
	tot_rama += resE[rama]
	tot_omega += resE[omega]
	tot_p_aa_pp += resE[p_aa_pp]
	tot_hbond_sr_bb += resE[hbond_sr_bb]
	tot_hbond_lr_bb += resE[hbond_lr_bb]
	tot_hbond_bb_sc += resE[hbond_bb_sc]
	tot_hbond_sc += resE[hbond_sc]
	tot_pro_close += resE[pro_close]
	tot_ref += resE[ref]
	
#, resE[fa_rep], resE[fa_sol], resE[fa_intra_rep], resE[fa_elec], resE[fa_dun], resE[dslf_fa13], resE[rama], resE[omega], resE[p_aa_pp], resE[hbond_sr_bb], resE[hbond_lr_bb], resE[hbond_bb_sc], resE[hbond_sc], resE[pro_close], resE[ref]
	totE = resE[fa_atr] + resE[fa_rep] + resE[fa_sol] + resE[fa_intra_rep] + resE[fa_elec] + resE[fa_dun] + resE[dslf_fa13] + resE[rama] + resE[omega] + resE[p_aa_pp] + resE[hbond_sr_bb] + resE[hbond_lr_bb] + resE[hbond_bb_sc] + resE[hbond_sc] + resE[pro_close] # + resE[ref]
	print("%4d %11f %11f %11f %12f %11f %11f %11f %11f %11f %11f %11f %11f %11f %11f %11f %11f %11f") % (i, resE[fa_atr], resE[fa_rep], resE[fa_sol], resE[fa_intra_rep], resE[fa_elec], resE[fa_dun], resE[dslf_fa13], resE[rama], resE[omega], resE[p_aa_pp], resE[hbond_sr_bb], resE[hbond_lr_bb], resE[hbond_bb_sc], resE[hbond_sc], resE[pro_close], resE[ref], totE)
	res_energy += totE

# Weighting contributions
tot_fa_atr *= scorefxn[fa_atr]
tot_fa_rep *= scorefxn[fa_rep]
tot_fa_sol *= scorefxn[fa_sol]
tot_fa_intra_rep *= scorefxn[fa_intra_rep]
tot_fa_elec *= scorefxn[fa_elec]
tot_fa_dun *= scorefxn[fa_dun]
tot_dslf_fa13 *= scorefxn[dslf_fa13]
tot_rama *= scorefxn[rama]
tot_omega *= scorefxn[omega]
tot_p_aa_pp *= scorefxn[p_aa_pp]
tot_hbond_sr_bb *= scorefxn[hbond_sr_bb]
tot_hbond_lr_bb *= scorefxn[hbond_lr_bb]
tot_hbond_bb_sc *= scorefxn[hbond_bb_sc]
tot_hbond_sc *= scorefxn[hbond_sc]
tot_pro_close *= scorefxn[pro_close]
tot_ref *= scorefxn[ref]

tot_resE = tot_fa_atr + tot_fa_rep + tot_fa_sol + tot_fa_intra_rep + tot_fa_elec + tot_fa_dun + tot_dslf_fa13 + tot_rama + tot_omega + tot_p_aa_pp + tot_hbond_sr_bb + tot_hbond_lr_bb + tot_hbond_bb_sc + tot_hbond_sc + tot_pro_close # + tot_ref

print("Total residue energy unweighted, res_energy= %f\n") % res_energy

print("tot_fa_atr= %11f") % tot_fa_atr
print("tot_fa_rep = %11f") % tot_fa_rep
print("tot_fa_sol = %11f") % tot_fa_sol
print("tot_fa_intra_rep = %11f") % tot_fa_intra_rep
print("tot_fa_elec = %11f") % tot_fa_elec
print("tot_fa_dun = %11f") % tot_fa_dun
print("tot_dslf_fa13 = %11f") % tot_dslf_fa13
print("tot_rama = %11f") % tot_rama
print("tot_omega = %11f") % tot_omega
print("tot_p_aa_pp = %11f") % tot_p_aa_pp
print("tot_hbond_sr_bb = %11f") % tot_hbond_sr_bb
print("tot_hbond_lr_bb = %11f") % tot_hbond_lr_bb
print("tot_hbond_bb_sc = %11f") % tot_hbond_bb_sc
print("tot_hbond_sc = %11f") % tot_hbond_sc
print("tot_pro_close = %11f") % tot_pro_close
print("tot_ref = %11f") % tot_ref

print("\ntot_resE = %11f") % tot_resE

print(scorefxn.show(pose))

# for i in range( rri , rrf+1 ):
#	print pose.energies().show(i)

end = time.time()
print("Energy time: %f") % (end - start)

#packer_task(pose, PDB_out)

# pose2 = rosetta.protocols.grafting.return_region(pose,rri-1,rrf+1)
# pose2.dump_pdb('myloop.pdb')


