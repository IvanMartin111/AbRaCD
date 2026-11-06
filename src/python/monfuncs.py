#!/usr/bin/env python

##!/usr/bin/python3

#import rosetta.core.pack.task    # for using resfiles
#from rosetta import *

import time
import pyrosetta
import pyrosetta.rosetta as rosetta
from pyrosetta import Pose, standard_packer_task, standard_task_factory, Vector1, MoveMap #, pose_from_file, ScoreFunction, create_score_function, get_fa_scorefxn
from pyrosetta.rosetta import core, protocols

#import rosetta.utility   # for using resfiles
#from rosetta.utility import *


################################################################################
# A GENERAL EXPLANATION
#!/home/mon/prgs/PyRosetta.Ubuntu-12.04LTS.64Bit.monolith/ipython.py
#!usr/bin/env python


# The AutoVivification function is to easily deal with hashes-of-hashes
# http://stackoverflow.com/questions/651794/whats-the-best-way-to-initialize-a-dict-of-dicts-in-python
class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

"""
Bla, bla, bla...
"""
def getKey(item):
	# print item[0]
	return float(item[0])

def getKey2(item):
	# print item[3]
	return float(item[3])

def keepbest(nbest,pose_list,pose,score,data):

	if len(pose_list) < nbest:
		pose_list.append([float(score),pose,data,0,0]) # Add score/pose/data
	else:
		for i in range(0,len(pose_list)):
			if pose_list[i][0] > score:  # substitute some worst element element, if necessary
				pose_list[i][0] = float(score)
				pose_list[i][1] = pose
				pose_list[i][2] = data
				return  # stops as soon as any worst element is found!

def loop_packer(test_pose,scorefxn,resfile):
    """
    Bla, bla...
    """
    # create a copy of the pose
#    test_pose = Pose()
#    test_pose.assign(pose)

    # this object is contained in PyRosetta v2.0 and above
#    pymover = PyMOL_Mover()

    # since the PackerTask specifies how the sidechains change, it has been
    #    extended to include sidechain constitutional changes allowing
    #    protein design, this method of design is very similar to sidechain
    #    packing; all rotamers of the possible mutants at a single residue
    #    are considered and the lowest scoring conformation is selected
    # design options include:
    #    -allow all amino acids
    #    -allow all amino acids except cysteine
    #    -allow specific amino acids
    #    -prevent specific amino acids
    #    -allow polar amino acids only
    #    -prevent polar amino acids
    #    -allow only the native amino acid
    # the myriad of packing and design options can be set manually or, more
    #    commonly, using a specific file format known as a resfile
    #    resfile syntax is explained at:
    #    http://www.rosettacommons.org/manuals/archive/rosetta3.1_user_guide/file_resfiles.html
    # manually setting deign options is tedious, the methods below are handy
    #    for creating resfiles

    # setup the design PackerTask, use the generated resfile
    pose_design = standard_packer_task(test_pose)
   # pose_design.or_optimize_h_mode(True) # if "True", turns on optimize_H_mode for all residues
    rosetta.core.pack.task.parse_resfile(test_pose, pose_design, resfile )
 #   pose_design.or_optimize_h_mode(True) # if "True", turns on optimize_H_mode for all residues
 #   print pose_design

    # prepare a new structure
    #test_pose.assign(pose)

    # perform design
    designmover = protocols.minimization_packing.PackRotamersMover(scorefxn, pose_design)
#    designmover = MinPackMover(scorefxn, pose_design) # similar to PackRotamersMover but with additional "off-rotamer" minimization steps (100 times slower)

 #   print 'Pre-design score:', scorefxn(test_pose)
 #   print 'Pre-design sequence: ...' + test_pose.sequence()[ri-1:rf] + '...'

    designmover.apply(test_pose)    # perform design
#    print '\nPost-design score:', scorefxn(test_pose)
#    print 'Post-design sequence: ...' + \
#        test_pose.sequence()[ri-1:rf] + '...'
#    test_pose.pdb_info().name( 'designed' )    # for PyMOL_Mover
#    pymover.apply(test_pose)

# NEIGHBORHOOD detection routine
#  pose --> Pose to repack
# Optional: 
#  dist --> Distance cutoff (<0: Full-pose repacking, =0: Loop-only repacking not-including anchors, >0: Distance dependant repacking)
#  rri ---> N-t anchor residue index (Rosetta numbering schemes for residue, use info.pdb2pose function)
#  rrf ---> C-t anchor residue index (Rosetta numbering schemes for residue, use info.pdb2pose function)
#  neighs --> Hash with the neighbor indices (if !=0 the neighborhood will be incremental)
# OUTPUT:  hash of residue indices in Rosetta numeration scheme
def neighbors(pose,dist=0,rri=-1,rrf=-1,neighs=0):
	t0 = time.time()
	if neighs == 0:
		neighs = {} # Hash (dictionary) with the neighbors indices

	if dist < 0 or rri < 0 or rrf < 0:
		print("neighbors()> Full protein repacking...")
		for r in range( 1, pose.total_residue()+1 ):  # residue indices in Rosetta numeration scheme (not-including anchors)
			res = pose.residue(r)
			if "CYD" not in res.name():
				# print "Residue %d is not CYD --> selected for repacking" % (r)
				#task_pack.temporarily_set_pack_residue(r, True)
				neighs[r]=1  # append residue index

	else:
		if dist == 0:  # Standard "just-loop" repacking (not-including anchors)

			print("neighbors()> Just-loop repacking...")

			# Screen loop residues
			#for r in range( rri , rrf+1 ):  # residue indices in Rosetta numeration scheme
			for r in range( rri+1 , rrf ):  # residue indices in Rosetta numeration scheme (not-including anchors)
				res = pose.residue(r)
				if "CYD" not in res.name():
					# print "Residue %d is not CYD --> selected for repacking" % (r)
					#task_pack.temporarily_set_pack_residue(r, True)
					neighs[r]=1  # append residue index
				
		else:  # Loop and neighborhood repacking
	
			print("neighbors()> Loop and %f A neighborhood repacking..." % (dist))

			dist2 = dist * dist  # Squared distance is faster...

			# Screen loop residues
			for r in range( rri+1 , rrf ):  # residue indices in Rosetta numeration scheme (without anchors)
				if "GLY" not in pose.residue(r).name():
					atype='CB'  # Non-GLY aminoacids have CB
				else:
					atype='CA'  # GLY only has CA instead
				atom = pose.residue(r).atom(atype)  # get CA/CB atom "r"
				xyz1 = atom.xyz()
				for x in range(1,pose.total_residue()+1): # screen all residues
					res = pose.residue(x)
					# print "AA: %s" % (res.name())
					if "GLY" not in res.name():
						atom2 = res.atom('CB')  # get CB atom "x"
						xyz2 = atom2.xyz()
						d2 = (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 + (xyz1[2]-xyz2[2])**2
						if d2 < dist2 and "CYD" not in res.name():
							# print "Residues %d and %d are neigbors (d2=%f) and are not CYD" % (r,x,d2)
							#task_pack.temporarily_set_pack_residue(x, True)
							neighs[x]=1  # append residue index
					else:
						atom2 = res.atom('CA')  # get CA atom "x" (Mon 25/1/2016 fixed: GLYs were not selected before!)
						xyz2 = atom2.xyz()
						d2 = (xyz1[0]-xyz2[0])**2 + (xyz1[1]-xyz2[1])**2 + (xyz1[2]-xyz2[2])**2
						if d2 < dist2:
							# print "Residues %d and %d are neigbors (d2=%f) and are not CYD" % (r,x,d2)
							#task_pack.temporarily_set_pack_residue(x, True)
							neighs[x]=1  # append residue index

	print("Total residues: %d  Neighboring residues: %d  (loop included)" % (pose.total_residue(),len(neighs)))
	tf = time.time()
	print("Neighborhood detection time: %f s" % (tf-t0))
	return neighs  # Output neighbors

#print pose.pdb_info().chain(500)
#print pose.pdb_info().number(500)

# Translate selection mask from Rosetta indices to PDB indices
#  neighs --> hash of residue indices in Rosetta numeration scheme
#  pose --> Pose to get PDB indices from
# OUTPUT:  hash of residue indices and chain ids in PDB numeration scheme
def rosseta2pdb(neighs,pose):
	t0 = time.time()
#	neighs2 = {} # Hash (dictionary) with the neighbors indices in PDB numeration
	neighs2 = AutoVivification() # Hash (dictionary) with the neighbors indices in PDB numeration
	info = pose.pdb_info()

	#Lookup the Rosetta internal number for residue 100 of chain A:
	#print pose.pdb_info().pdb2pose('A', 100)
	# The converse command is:
	#print pose.pdb_info().pose2pdb(25)
	for x in neighs:  # residue indices in Rosetta numeration scheme
		index, chain = info.pose2pdb(x).split()
#		print "%d %s %s chain=%s %d" % (x,info.pose2pdb(x).split(),index,chain, int(index))
#		neighs2[ int(index) ] = chains  # append residue index
		neighs2[ int(index) ][ chain ] = 1  # append residue index

	return neighs2

# Translate selection mask from PDB indices to Rosseta indices
#  neighs --> hash of residue indices in PDB numeration scheme
#  pose --> Pose to get Rosseta indices from
#  chain --> Chain-ID (one-letter) required
# OUTPUT:  hash of residue indices in Rosetta numeration scheme
def pdb2rosseta(neighs,pose,chain):
	t0 = time.time()
	neighs2 = {} # Hash (dictionary) with the neighbors indices in PDB numeration
	info = pose.pdb_info()

	#Lookup the Rosetta internal number for residue 100 of chain A:
	#print pose.pdb_info().pdb2pose('A', 100)
	# The converse command is:
	#print pose.pdb_info().pose2pdb(25)
	for x in neighs:  # residue indices in Rosetta numeration scheme (key values)
		#print "myindex= %d" % int(info.pdb2pose(chain,x))
		for c in neighs[x]: # Chains with residue index "x"
			index = int(info.pdb2pose(c,x))
			if index != 0:
				neighs2[ index ] = 1  # append residue index
			else:
				print("Warning: residue %d not found!" % x)

	return neighs2


# Translate selection mask from PDB indices to Rosseta indices
#  neighs --> hash of residue indices in PDB numeration scheme
#  pose --> Pose to get Rosseta indices from
#  chain --> Chain-ID (one-letter) required
# OUTPUT:  hash of residue indices in Rosetta numeration scheme
def pdb2rosseta_old(neighs,pose,chain):
	t0 = time.time()
	neighs2 = {} # Hash (dictionary) with the neighbors indices in PDB numeration
	info = pose.pdb_info()

	#Lookup the Rosetta internal number for residue 100 of chain A:
	#print pose.pdb_info().pdb2pose('A', 100)
	# The converse command is:
	#print pose.pdb_info().pose2pdb(25)
	for x in neighs:  # residue indices in Rosetta numeration scheme
		#print "myindex= %d" % int(info.pdb2pose(chain,x))
		if int(info.pdb2pose(chain,x)) != 0:
			neighs2[ int(info.pdb2pose(chain,x)) ]=1  # append residue index
		else:
			print("Warning: residue %d not found!" % x)

	return neighs2


# Loop and Neighborhood repacking routine
#  pose --> Pose to repack
#  scorefxn --> Energy function
#  myneighbors --> Neighborhood mask (Rosseta numbering)
#  task_pack --> Repacking object (it should be already initialized)
#  dist --> Distance cutoff
#def neigbor_packer(pose,scorefxn,task_pack,dist):
#
# Good baseline settings for side chain packing:  https://www.rosettacommons.org/node/3983
# Repacking options: https://www.rosettacommons.org/manuals/rosetta3.1_user_guide/opt_packing.html
def neighbor_packer(pose,scorefxn,myneighbors):

	t0 = time.time()
	# Repacking object stuff...
	print("About to create task_pack...")
	task_pack = standard_packer_task(pose)  # create object
	# By default, only rotamers are considered for packing and the original sidechain conformation will be lost.
	# See "Packing" documentation: http://www.pyrosetta.org/documentation
	print("About to create restrict_to_repacking...")
	task_pack.restrict_to_repacking()  # only repacking
#	task_pack.temporarily_fix_everything()  # fix all first

#	info = pose.pdb_info()


#	for x in range(1,pose.total_residue()+1): # screen all residues
#		res = pose.residue(x)
#		# print "AA: %s" % (res.name())
#		if "CYS" not in res.name():
#			task_pack.temporarily_set_pack_residue(x, False)


#	vector1 = rosetta.utility.vector1_boolean()
#	vector1 = rosetta.Vector1([])

	v = []
	for x in range(1,pose.total_residue()+1): # screen all residues
#		vector1.append(True)
		v.append(False)


#	for i in vector1:
#		print "mierda %d ," % i

#	print vector1

	for x in myneighbors:  # residue indices in Rosetta numeration scheme
		#print "x=%d --> %s" % (x,info.pose2pdb(x))
		v[x-1] = True
		# This ("old" and "temporary") function seems broken when accessing CYS..., see:
		# https://www.rosettacommons.org/content/individual-residue-repacking-causes-segmentation-fault-or-error-seqpos-1-error-exit-srccorec
		# task_pack.temporarily_set_pack_residue(x, True) 

	# print "Repacking mask: "
	# print v

	# This seems to be the way to create a PyRosetta's "vector1"... 
	# See PyRosetta's FAQ: http://www.pyrosetta.org/faq#TOC-2.-How-do-I-construct-Vector1-objects-
	# Examples:
	# print rosetta.Vector1( [ 1.0 , 2.0 , 3.0 ] )
	# print rosetta.Vector1( [ True , False , True ] )
	# print rosetta.Vector1( [ 'a' , 'b' , 'c' ] )
	vector1 = Vector1(v) 

	# print "My vector1"
	# print vector1
	task_pack.restrict_to_residues(vector1) # Apply repacking only to "True" residues, "False" will not be repacked

	# print task_pack


	# Repacking
	t0 = time.time()
	packmover = protocols.minimization_packing.PackRotamersMover(scorefxn,task_pack)
	packmover.apply(pose)
	tf = time.time()
	print("Repacking time: %f s" % (tf-t0))


def select_chi(movemap,selection):
	#movemap.set_chi_true_range(rri+1,rrf-1)
	for i in selection:
		movemap.set_chi(i,True)
#	else:
#		movemap.set_chi(i,False)


def write_res(f,pose,selection,bf):
	# for each residue
	ati = 1;
	info = pose.pdb_info()
	chain = info.chain(pose.residue(selection[0]).chain())
	for i in selection:
		resi = pose.residue(i) # get residue
		res_name = resi.name()
		if res_name.find('HIS_D') != -1: # PyRosetta evaluates the protonation state of HIS
			res_name = 'HIS'
		if res_name.find('CYD') != -1: # PyRosetta evaluates the oxidation state of Disulfide bonds (CYS or CYD)
			res_name = 'CYS'
		res_num = info.number(i)
		# for each atom
		for atom in range( resi.natoms() ):
			atom_name = resi.atom_name( atom + 1 )
			if atom_name.find('OV') != -1: # PyRosetta introduces "OV" duplicate atoms during ChainTree creation
				continue
			coords = resi.xyz( atom + 1 )
			f.write('ATOM  %5d %-5s%3s %c%4d    %8.3lf%8.3lf%8.3lf%6.2f%6.2f\n' % (ati, atom_name, res_name[0:3], chain, res_num, coords[0], coords[1], coords[2], 1.0, bf[i-1] ) )
			ati += 1 # atom index

def write_pdb(pose, myfile, bf):
	f=open(myfile, 'w+')
	selection = range( 1 , pose.total_residue() + 1 )
	write_res(f,pose,selection,bf) # Write all residues in selection with PDB-format
	f.close()

def append_mpdb(pose, myfile, bf, model):
	f=open(myfile, 'a')
	f.write('MODEL %6d\n' % model)
	selection = range( 1 , pose.total_residue() + 1 )
	write_res(f,pose,selection,bf) # Write all residues in selection with PDB-format
	f.write('ENDMDL\n')
	f.close()



# READ MULTI-PDB (untested yet)
#  outputs: pdbs --> Cartesian coordinates 4D array: pdbs[loop][atom][residue][coord]
def read_multipdb(pdb_loop):

	#pose, resfilename, pack = True, design = False, input_sc = True, freeze = [], specific = {}):
	# GET LOOP COORDINATES FROM MULTI-PDB
	pdbs = [] # 4D array with all loop coordinates
	coords_N = []
	coords_CA = []
	coords_C = []
	npdb = 0
	for line in open(pdb_loop):
		list = line.split()
		id = list[0]

		if id == 'ATOM':
			type = list[2]
			pos = [float(list[6]), float(list[7]), float(list[8])]
			# print pos
			if type == 'N':
			    coords_N.append(pos)
			if type == 'CA':
			    coords_CA.append(pos)
			if type == 'C':
			    coords_C.append(pos)
		elif id == 'ENDMDL':
			# Checking input
			if ( len(coords_N) != len(coords_CA) or len(coords_N) != len(coords_C) ):
			     print("Error, missing N, CA or C atoms...")
			     exit()
			if ( len(coords_N) != loop_len+2 ):
			     print("Error, mismatch between loop coordinates and residue indices")
			     exit()
			coords = []
			coords.append(coords_N)
			coords.append(coords_CA)
			coords.append(coords_C)
			pdbs.append(coords)
			coords_N = []
			coords_CA = []
			coords_C = []
			npdb += 1

	for l in range(0,len(pdbs)):
		print("Loop coordinates for N:")
		print(pdbs[l][0])
		print("Loop coordinates for CA:") 
		print(pdbs[l][1])
		print("Loop coordinates for C:")
		print(pdbs[l][2])

	print("Loaded %d loops from Multi-PDB %s" % (npdb, pdb_loop))
	return pdbs

# UPDATE BACKBONE (N, CA, C) COORDINATES and REPACK SIDECHAINS
#  pose --> PyRosetta's pose
#  pdbs --> Cartesian coordinates 4D array: pdbs[loop][atom][residue][coord]
def update_and_repack(pose,pdbs):
	for l in range(0,len(pdbs)):
		i = 1
		for r in range( ri , rf+1 ):
		     atom = pose.residue( r ).atom('N')
		     coord = atom.xyz()
		     coord.x = pdbs[l][0][i][0]
		     coord.y = pdbs[l][0][i][1]
		     coord.z = pdbs[l][0][i][2]
		     atom = pose.residue( r ).atom('CA')
		     coord = atom.xyz()
		     coord.x = pdbs[l][1][i][0]
		     coord.y = pdbs[l][1][i][1]
		     coord.z = pdbs[l][1][i][2]
		     atom = pose.residue( r ).atom('C')
		     coord = atom.xyz()
		     coord.x = pdbs[l][2][i][0]
		     coord.y = pdbs[l][2][i][1]
		     coord.z = pdbs[l][2][i][2]
		     i += 1 # residue counter

		start = time.time()
		loop_packer(pose,scorefxn)
		end = time.time()
		print("Loop_packer time: %f" % (end - start))
		energy = scorefxn(pose)
		print("Loop %5d energy: %f" % (l, energy))

		# Dump loop into Pose() object
		loop = rosetta.protocols.grafting.return_region(pose,ri-1,rf+1)
		print("Storing loop %d, residues from %d to %d" % (l+1, ri-1, rf+1))
		# Apply PDBinfo from dummy...
		loop.pdb_info(loop_info)
		loop_name = "myloop_%04d.pdb" % (l+1)
		loop.dump_pdb(loop_name)

# READ DIHEDRALS from table with this format:
# Omega_1, Phi_1, Psi_1, Omega_2, ..., Omega_N, Phi_N, Psi_N
def read_dihedrals(dh_file):
	debug = 0
	n_new  = 0
	n_old = 0
	index = 0
	dihedrals = [] # 1D array with the dihedral angles of all loops

	for line in open(dh_file):
		list = line.split()
		if(list[0] != '#'):
			n_old = n_new
			n_new = len(list)
			if(index != 0):
				if(n_new != n_old):
					print("Error, number of dihedrals mismatch (n_new != n_old, %d != %d) in file: %s\n" % (n_new, n_old, dh_file))
					exit()
			floats = [float(x) for x in list] # convert string to float
			dihedrals.append(floats) # Dihedrals array: Omega_1, Phi_1, Psi_1, Omega_2, ..., Omega_N, Phi_N, Psi_N
			index += 1 # Counts the number of lines, i.e. the number of loops
	if debug:		
		for l in range(0,len(dihedrals)):
			print("Dihedral angles for loop %d" % (l+1))
			print(dihedrals[l])

	print("Readed %d dihedral angles from %d loops (form %s)" % (n_new, index, dh_file))
	return dihedrals

# READ VALENCE ANGLES from table with this format:
# NCAC_1, CACN_1, CNCA_1, NCAC_2, ..., NCAC_3, CACN_3, CNCA_3, NCAC_N. 
# Note that the final CACN_N and CNCA_N angles should not be provided.
def read_valences(val_file):
	debug = 0
	n_new  = 0
	n_old = 0
	index = 0
	valences = [] # 1D array with the valence angles of all loops

	for line in open(val_file):
		list = line.split()
		if(list[0] != '#'):
			n_old = n_new
			n_new = len(list)
			if(index != 0):
				if(n_new != n_old):
					print("Error, number of valence angles mismatch (n_new != n_old, %d != %d) in file: %s\n" % (n_new, n_old, val_file))
					exit()
			floats = [float(x) for x in list] # convert string to float
			valences.append(floats) # Valence angles array: NCAC_1, CACN_1, CNCA_1, NCAC_2, ..., NCAC_3, CACN_3, CNCA_3, NCAC_N.
			index += 1 # Counts the number of lines, i.e. the number of loops

	if debug:		
		for l in range(0,len(valences)):
			print("Valence angles for loop %d" % (l+1))
			print(valences[l])

	print("Readed %d valence angles from %d loops (form %s)" % (n_new, index, val_file))
	return valences

# READ BOND LENGTHS from table with this format:
# N-CA_1, CA-C_1, C-N_1, N-CA_2, ..., N-CA_3, CA-C_3, C-N_3, N-CA_4, CA-C_4.
# Note that the final C-N bond length should not be provided.
def read_lengths(len_file):
	debug = 0
	n_new  = 0
	n_old = 0
	index = 0
	lengths = [] # 1D array with the bond lengths of all loops

	for line in open(len_file):
		list = line.split()
		if(list[0] != '#'):
			n_old = n_new
			n_new = len(list)
			if(index != 0):
				if(n_new != n_old):
					print("Error, number of bond lengths mismatch (n_new != n_old, %d != %d) in file: %s\n" % (n_new, n_old, len_file))
					exit()
			floats = [float(x) for x in list] # convert string to float
			lengths.append(floats) # Bond lengths array: N-CA_1, CA-C_1, C-N_1, N-CA_2, ..., N-CA_3, CA-C_3, C-N_3, N-CA_4, CA-C_4.
			index += 1 # Counts the number of lines, i.e. the number of loops

	if debug:		
		for l in range(0,len(lengths)):
			print("Bond lengths for loop %d" % (l+1))
			print(lengths[l])

	print("Readed %d bond lengths from %d loops (form %s)" % (n_new, index, len_file))
	return lengths

# Reads some column from file
#  rmsds_file --> Input plain-text file
#  col --> Column index (0,1,2,...)
#  ignore --> Lines that start with this word will be discarded
def read_column(rmsds_file,col,ignore):
	rmsds = [] # 1D array with the bond lengths of all loops
	for line in open(rmsds_file):
		list = line.split()
		if(list[0] != ignore):
			rmsds.append(list[col]);
	return rmsds

# Get dihedral angles from pose
def get_dihedrals(pose, resi=1, nres=None): # pyrosetta starts at 1
	dhs = []
	if nres is None:
		nres = pose.total_residue()
	for l in range(0,nres): # iters from resi to resi+nres residues
		dhs.append(pose.phi(resi+l)) # Phi angle
		dhs.append(pose.psi(resi+l)) # Psi angle
	return dhs

def update_dihedrals(pose, dhs, resi):
	rad2deg = float(180/3.14159265358979);
	nres = int(len(dhs)/3); # number of loop residues
	for l in range(0,nres): # iters from resi to resi+nres residues
		# print "Updating dihedrals for residue %d (nres= %d)" % (resi+l,nres)
		# Omega_1, Phi_1, Psi_1, Omega_2, ..., Omega_N, Phi_N, Psi_N
		pose.set_phi(resi+l,dhs[l*3+1]*rad2deg) # Phi angle
		pose.set_psi(resi+l,dhs[l*3+2]*rad2deg) # Psi
		if (l < nres-1): # Non-last omega
			pose.set_omega(resi+l,dhs[l*3+3]*rad2deg) # check this

def update_valences(pose, vals, resi):
	rad2deg = float(180/3.14159265358979);
	nres = int((len(vals)+2)/3); # number of loop residues

	at_N = core.id.AtomID(1,resi)
	at_CA = core.id.AtomID(2,resi)
	at_C = core.id.AtomID(3,resi)
	pose.conformation().set_bond_angle(at_N,at_CA,at_C,vals[0]) # setting NCAC_1
	
	for l in range(1,nres): # iters from resi+1 to resi+nres residues
		conf = pose.conformation()
		# print "Updating valence angles for residue %d (nres= %d)" % (resi+l,nres)
		# NCAC_1, CACN_1, CNCA_1, NCAC_2, ..., NCAC_3, CACN_3, CNCA_3, NCAC_N.
		at_N1 = core.id.AtomID(1,resi+l)
		at_CA1 = core.id.AtomID(2,resi+l)
		at_C1 = core.id.AtomID(3,resi+l)
		conf.set_bond_angle(at_CA,at_C,at_N1,vals[(l-1)*3+1]) # setting CACN_(N-1)
		conf.set_bond_angle(at_C,at_N1,at_CA1,vals[(l-1)*3+2]) # setting CNCA_(N-1)
		conf.set_bond_angle(at_N1,at_CA1,at_C1,vals[l*3]) # setting NCAC_N
		at_N = at_N1 # store previous N
		at_CA = at_CA1 # store previous CA
		at_C = at_C1 # store previous C

def update_lengths(pose, lens, resi):
	nres = int((len(lens)+1)/3); # number of loop residues

	at_N = core.id.AtomID(1,resi)
	at_CA = core.id.AtomID(2,resi)
	at_C = core.id.AtomID(3,resi)
	pose.conformation().set_bond_length(at_N,at_CA,lens[0]) # setting N-CA_1
	pose.conformation().set_bond_length(at_CA,at_C,lens[1]) # setting CA-C_1

	for l in range(1,nres): # iters from resi+1 to resi+nres residues
		conf = pose.conformation()
		# print "Updating bond lengths for residue %d (nres= %d)" % (resi+l,nres)
		# N-CA_1, CA-C_1, C-N_1, N-CA_2, ..., N-CA_3, CA-C_3, C-N_3, N-CA_4, CA-C_4.
		at_N1 = core.id.AtomID(1,resi+l)
		at_CA1 = core.id.AtomID(2,resi+l)
		at_C1 = core.id.AtomID(3,resi+l)
		conf.set_bond_length(at_C,at_N1,lens[l*3-1]) # setting C-N_(N-1)
		conf.set_bond_length(at_N1,at_CA1,lens[l*3]) # setting N-CA_(N)
		conf.set_bond_length(at_CA1,at_C1,lens[l*3+1]) # setting CA-C_(N)
		at_N = at_N1 # store previous N
		at_CA = at_CA1 # store previous CA
		at_C = at_C1 # store previous C

# Compute RMSD between two poses
#  pose1 --> First pose
#  pose2 --> Second pose (same size as first pose)
#  select --> residue indices to be considered (selection)
#  consider --> Atom IDs: 'C', 'N', 'CA', etc... (Python's list of strings)
#  discard --> Atoms IDs matching substrings in "discard" list will be neglected. (To discard Hydrogens: 'H')
def rmsd(pose1,pose2,select=None,consider=None,discard=None):
	dx = 0.0;
	dy = 0.0;
	dz = 0.0;
	natom = 0; # number of atoms counter

	if select is None:
		select = range(1,pose1.total_residue()) # then select all residues

	if discard is None and consider is not None: # When no-atoms shoud be discarded and consider exists
		for i in select:
			resi = pose1.residue(i) # get residue 1
			for atom in range( resi.natoms() ): # iter residue atoms
				atom_name = resi.atom_name( atom + 1 ).strip()
				for con in consider:
					if atom_name == con:
						resj = pose2.residue(i) # get residue 2
						for atom2 in range( resj.natoms() ): # iter residue atoms
							atom_name2 = resj.atom_name( atom2 + 1 ).strip()
							#print("Atom %s %s %s \n" % (i,atom_name,atom_name2))
							if atom_name2 == con:
								xyz1 = resi.xyz( atom + 1 ) # get coords 1
								xyz2 = resj.xyz( atom2 + 1 ) # get coords 2
								dx += (xyz1[0]-xyz2[0])**2 
								dy += (xyz1[1]-xyz2[1])**2 
								dz += (xyz1[2]-xyz2[2])**2

								natom += 1 # counting the number of atoms
								
								#print("Atom %s %s %f %f %f  -> %f %f %f\n" % (atom_name,atom_name2, xyz1[0],xyz1[1],xyz1[2],xyz2[0],xyz2[1],xyz2[2]))

								break # exit for loop after atom is found
	elif consider is None and discard is not None: # when some atoms should be discarded only
		ndis = len(discard)
		for i in select:
			resi = pose1.residue(i) # get residue 1
			for atom in range( resi.natoms() ): # iter residue atoms
				atom_name = resi.atom_name( atom + 1 ).strip()
				neglect = False
				for dis in discard:
					# if atom_name == neg
					#print("Discard  %s %s %s\n"  %(dis,atom, atom_name ))
					if atom_name.find(dis,0,len(atom_name)) != -1:
						neglect = True
						break
				if not neglect:
					#print("ACEPT  %s %s\n"  %(atom, atom_name ))
					resj = pose2.residue(i) # get residue 2
					for atom3 in range( resj.natoms() ): # iter residue atoms
						atom_name3 = resj.atom_name( atom3 + 1 ).strip()
						#print("Atom %d -%s-%s-" % (i,atom_name,atom_name3))
						if atom_name3 == atom_name:
							xyz1 = resi.xyz( atom + 1 ) # get coords 1
							xyz2 = resj.xyz( atom3 + 1 ) # get coords 2
							dx += (xyz1[0]-xyz2[0])**2 
							dy += (xyz1[1]-xyz2[1])**2 
							dz += (xyz1[2]-xyz2[2])**2 
							natom += 1 # counting the number of atoms
							#print("Atom %s %s %f %f %f  -> %f %f %f" % (atom_name,atom_name3, xyz1[0],xyz1[1],xyz1[2],xyz2[0],xyz2[1],xyz2[2]))
							break # exit for loop after atom is found
	elif consider is not None and discard is not None: # when atoms should be discarded and considered
		ndis = len(discard)
		for i in select:
			resi = pose1.residue(i) # get residue 1
			for atom in range( resi.natoms() ): # iter residue atoms
				atom_name = resi.atom_name( atom + 1 ).strip()
				atom_name2 = atom_name.strip()
				for con in consider:
					neglect = False
					for dis in discard:
						# if atom_name == neg
						if atom_name.find(dis,0,len(atom_name)) != -1:
							neglect = True
							break
					if atom_name2 == con and not neglect:
						resj = pose2.residue(i) # get residue 2
						for atom3 in range( resj.natoms() ): # iter residue atoms
							atom_name3 = resj.atom_name( atom3 + 1 ).strip()
							#print("Atom %s %s %s \n" % (i,atom_name,atom_name2))
							if atom_name3 == atom_name:
								xyz1 = resi.xyz( atom + 1 ) # get coords 1
								xyz2 = resj.xyz( atom3 + 1 ) # get coords 2
								dx += (xyz1[0]-xyz2[0])**2 
								dy += (xyz1[1]-xyz2[1])**2 
								dz += (xyz1[2]-xyz2[2])**2 
								natom += 1 # counting the number of atoms
								print("Atom %s %s %f %f %f  -> %f %f %f\n" % (atom_name,atom_name2, xyz1[0],xyz1[1],xyz1[2],xyz2[0],xyz2[1],xyz2[2]))
								break # exit for loop after atom is found
	else: # when all atoms should be considered
		for i in select:
			resi = pose1.residue(i) # get residue 1
			for atom in range( resi.natoms() ): # iter residue atoms
				atom_name = resi.atom_name( atom + 1 ).strip()
				resj = pose2.residue(i) # get residue 2
				for atom2 in range( resj.natoms() ): # iter residue atoms
					atom_name2 = resj.atom_name( atom2 + 1 ).strip()
					#print("Atom %s %s %s \n" % (i,atom_name,atom_name2))
					if atom_name2 == atom_name:
						xyz1 = resi.xyz( atom + 1 ) # get coords 1
						xyz2 = resj.xyz( atom2 + 1 ) # get coords 2
						dx += (xyz1[0]-xyz2[0])**2 
						dy += (xyz1[1]-xyz2[1])**2 
						dz += (xyz1[2]-xyz2[2])**2 
						natom += 1 # counting the number of atoms 
						#print("Atom %s %s %f %f %f  -> %f %f %f\n" % (atom_name,atom_name2, xyz1[0],xyz1[1],xyz1[2],xyz2[0],xyz2[1],xyz2[2]))
						break

	if natom == 0:
		return 0  # If for some reason there are not any atoms... (An all-GLY loop that has not sidechain atoms)
	else:
		return ( ( (dx + dy + dz)/natom )**0.5 ) # returning RMSD

#	return ( ( (dx + dy + dz)/natom )**0.5 ) # returning RMSD

# Loop relaxation script: corrected_relax_loops_in_situ.py
# From PD2 paper: MacDonald et al. PLOS (2013)
#  pose --> Input pose
def pd2relax(pose,scorefxn,fa_rep_wt,chbrk_wt,minmover,niters):

	scorefxn.set_weight(core.scoring.fa_rep,fa_rep_wt)
#        scorefxn.set_weight(core.scoring.chainbreak,chbrk_wt)
	scorefxn.set_weight(core.scoring.chainbreak,1.0)
	best_score = scorefxn(pose)
	best_pose = Pose()
	best_pose.assign(pose)

	for x in range(1, niters+1):
		#scorefxn.set_weight(core.scoring.fa_rep,0.02*fa_rep_wt)
		scorefxn.set_weight(core.scoring.fa_rep,0.01*fa_rep_wt)
		#scorefxn.set_weight(core.scoring.chainbreak,100.0)
		#scorefxn.set_weight(core.scoring.chainbreak,1.0)
		t0 = time.time()
		loop_packer(pose,scorefxn,resfile_name)
		tf = time.time()
		print("1st packer time: %f" % (tf-t0))
		score1 = scorefxn(pose)
		minmover.tolerance(0.01)
		minmover.score_function(scorefxn)
		minmover.apply(pose)
		score2 = scorefxn(pose)
		print("Repuslive weight %.2f. Energy from %10.3f to %10.3f" % (0.01, score1, score2))

		#scorefxn.set_weight(core.scoring.fa_rep,0.250*fa_rep_wt)
		scorefxn.set_weight(core.scoring.fa_rep,0.3*fa_rep_wt)
		#scorefxn.set_weight(core.scoring.chainbreak,100.0)
		#scorefxn.set_weight(core.scoring.chainbreak,1.0)
		loop_packer(pose,scorefxn,resfile_name)
		score1 = scorefxn(pose)
		minmover.tolerance(0.01)
		minmover.score_function(scorefxn)
		minmover.apply(pose)
		score2 = scorefxn(pose)
		print("Repuslive weight %.2f. Energy from %10.3f to %10.3f" % (0.3, score1, score2))
	
		#scorefxn.set_weight(core.scoring.fa_rep,0.550*fa_rep_wt)
		scorefxn.set_weight(core.scoring.fa_rep,0.6*fa_rep_wt)
		#scorefxn.set_weight(core.scoring.chainbreak,100.0)
		#scorefxn.set_weight(core.scoring.chainbreak,1.0)
		loop_packer(pose,scorefxn,resfile_name)
		score1 = scorefxn(pose)
		#minmover.tolerance(0.01)
		minmover.tolerance(0.001)
		minmover.score_function(scorefxn)
		minmover.apply(pose)
		score2 = scorefxn(pose)
		print("Repuslive weight %.2f. Energy from %10.3f to %10.3f" % (0.6, score1, score2))
	
		scorefxn.set_weight(core.scoring.fa_rep,fa_rep_wt)
		# scorefxn.set_weight(core.scoring.chainbreak,200.0)
		#scorefxn.set_weight(core.scoring.chainbreak,100.0)
		#scorefxn.set_weight(core.scoring.chainbreak,1.0)
		loop_packer(pose,scorefxn,resfile_name)
		score1 = scorefxn(pose)
		minmover.tolerance(0.00001)
		minmover.score_function(scorefxn)
		minmover.apply(pose)
		score2 = scorefxn(pose)
		print("Repuslive weight %.2f. Energy from %10.3f to %10.3f" % (1.0, score1, score2))
		
		# scorefxn.set_weight(core.scoring.chainbreak,chbrk_wt)
		# this_score = scorefxn(pose)
		this_score = score2
		if this_score <= best_score:
			best_pose.assign(pose)
			best_score = this_score
			print("Good score: " + str(best_score))
		else:
			print("Bad score: " + str(this_score))
			pose.assign(best_pose)

	print("Best score: " + str(best_score))
	return(best_score)


# "monrelax" - Loop relaxation function based on the script: 
#  corrected_relax_loops_in_situ.py (FastRelax) from PD2's paper: MacDonald et al. PLOS (2013)
#
#  pose --> Input pose
#  scorefxn --> Score function (the initial "fa_rep" weight will be restored before return)
#  task_pack --> Repacking object
#  niters --> Number of "outher minimization loops" 
#  thr --> Neighborhood threshold for repacking [A] ( <0 --> minimization neighborhood list will be used for repacking)
#  minneighbors --> Neighborhood (NB) list for MinMover
#  rri ---> N-t anchor residue index (Rosetta numbering schemes for residue, use info.pdb2pose function)
#  rrf ---> C-t anchor residue index (Rosetta numbering schemes for residue, use info.pdb2pose function)
#  repackouther --> =1 to Repack in outher loop, =0 to Repack in inner loop (default)
#  irepack --> =0 to repack in every inner loop iteration (VdW ramp), or choose some iteration to repack (1,2,...,N)
#  ntpsi --> =True to minimize N-t Psi dihedral angle, =False not minimize N-t Psi (default)
#  mintype --> Flavor of descent minimization for PyRosetta's MinMover ('dfpmin_armijo_nonmonotone' by default)
#              You can select among any \"Min mover\" minimization method of Rosetta: dfpmin (exact line search), dfpmin_armijo (inexact line search), dfpmin_armijo_nonmonotone (less exact line search), dfpmin_strong_wolfe (More-Thuente line minimization). See: https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d7/df1/minimization_overview.html for reference.
#  premature --> If True, then forces premature exit from outher loop.
#
# Note: ("pd2relax4" is now "monrelax")
def monrelax(pose,scorefxn,niters,thr,minneighbors,fas,tols,rri=-1,rrf=-1,repackouther=0,irepack=0,ntpsi=False,mintype='dfpmin_armijo_nonmonotone',premature=False):

	renb = False; # "re-neighborhood" selection

	# Back up initial weight
	fa_rep_wt = scorefxn.get_weight(core.scoring.fa_rep)
	# Get initial score
	best_score = scorefxn(pose)
	# Copy initial pose
	best_pose = Pose()
	best_pose.assign(pose)

	time_packing = 0
	time_minimization = 0

	# Required to force minimizer (MoveMap) to minimice Chi's !!!
	if thr >= 0 and rri >= 0 and rrf >= 0:
		repackneighbors = neighbors(pose,thr,rri,rrf)  # Define neighborhood for current pose
	else:
		repackneighbors = minneighbors # minimization neighborhood list will be used for repacking

	# "Outher-loop"
	for x in range(1, niters+1): 
		# Iter "fa" weights

		######################
		# REPACKING - OUTHER #
		######################
		#
		# Repack-outher enabled and avoiding innecessary repack given loop has been already repacked outside this function...
		if repackouther == 1:
			if irepack >= 0: # "irepack" has to be 0 or greater to perform any repacking
				if irepack == 0 or irepack == x: # repack always or only for selected iteration (1,2,...,niters)
					print("Repacking loop (thr=%f) at Outher loop (Outher=%d)\n" % (thr,x))
					if thr >= 0 and rri >= 0 and rrf >= 0: # Adaptative neighborhood for loop (and eventually, its environment)
						repackneighbors = neighbors(pose,thr,rri)  # Define neighborhood for repacking (adaptative)
					else:
						repackneighbors = minneighbors  # Define neighborhood for repacking equal to MinMover's NB
					neighbor_packer(pose,scorefxn,repackneighbors)  # Repack current pose

		# "Inner-loop"
		for i in range(len(fas)): 
			print("%d,%d) Repuslive weight %.2f and tolerance %f:" % (x, i, fas[i], tols[i]))
			scorefxn.set_weight(core.scoring.fa_rep,fas[i]*fa_rep_wt)

			#####################
			# REPACKING - INNER #
			#####################
			t0 = time.time()
			#out_name = "pre_relax4_%05d.pdb" % (i+1)
			#pose.dump_pdb(out_name)

			if repackouther == 0: # Inner repack enabled and avoiding innecessary repack 
				if irepack >= 0: # "irepack" has to be 0 or greater to perform any repacking
					if irepack == 0 or irepack == i+1: # repack always or only for selected iteration (1,2,...,niters)
						print("Repacking loop (thr=%f) in Inner loop (Outher=%d Inner=%d)\n" % (thr,x,i+1))
						if thr >= 0 and rri >= 0 and rrf >= 0: # Adaptative neighborhood for loop (and eventually, its environment)
							repackneighbors = neighbors(pose,thr,rri,rrf)  # Define neighborhood for repacking (adaptative)
						else:
							repackneighbors = minneighbors  # Define neighborhood for repacking equal to MinMover's NB
						neighbor_packer(pose,scorefxn,repackneighbors)  # Repack current pose

			score1 = scorefxn(pose)  # Energy before minimization

			#out_name = "post_relax4_%05d.pdb" % (i+1)
			#pose.dump_pdb(out_name)

			#print "pack1 %f" % (score1)
			#neighbor_packer(pose,scorefxn,minneighbors)  # Repack current pose
			#print "pack2 %f" % (scorefxn(pose))
			#neighbor_packer(pose,scorefxn,minneighbors)  # Repack current pose
			#print "pack3 %f" % (scorefxn(pose))

			tf = time.time()
			print("%d,%d) Packer time: %f" % (x, i, tf-t0))
			print("%d,%d) Packer Energy: %f , and time: %f" % (x, i, score1, tf-t0))
			time_packing += tf-t0

			# exit()

			################
			# MINIMIZATION #
			################
			#
			# The MinMover carries out a gradient based minimization to find the nearest local minimum in the energy function, 
			# such as that used in one step of the Monte-Carlo plus Minimization algorithm of Li & Scheraga.
			#
			t0 = time.time() # Minimization timer

			# This creates a new Fold-tree (it must be generated before MoveMap and MinMover for some reason...)
			if rri > 0 and rrf > 0:
				loop = protocols.loops.Loop(rri,rrf)
				loop.auto_choose_cutpoint(pose)
				protocols.loops.set_single_loop_fold_tree(pose,loop)
				protocols.loops.add_cutpoint_variants(pose)

			# Define variables to be minimized in the "MoveMap"
			mm4060 = MoveMap()
			select_chi(mm4060,minneighbors)

			if rri > 0 and rrf > 0: # ONLY minimize Backbone Phi's and Psi's when some loop is selected
				mm4060.set_bb_true_range(rri+1,rrf-1)  # <-- THIS MOVES BETWEEN ANCHORS (NOT-MOVES ANCHORS)
				mm4060.set(core.id.TorsionID(rrf, core.id.BB, 1), True)  # <-- Move C-terminal-anchor PHI dihedral angle (seems a bit faster)
				if ntpsi:
					mm4060.set(core.id.TorsionID(rri, core.id.BB, 2), True)  # <-- Move N-terminal-anchor PSI dihedral angle
			#mm4060.show()

			minmover = protocols.minimization_packing.MinMover() # Create MinMover
			minmover.movemap(mm4060) # Assign MoveMap to MinMover
			minmover.score_function(scorefxn) # Assign scorefunction to MinMover
			minmover.min_type(mintype) # Set minimization-type 
			# Check the following links to learn more about minimization-types:
			# https://www.rosettacommons.org/node/3395
			# Flavors of minimization in Rosetta:
			# https://www.rosettacommons.org/manuals/archive/rosetta3.5_user_guide/d7/df1/minimization_overview.html
			# "Jim recommends dfpmin_strong_wolfe or dfpmin_armijo_nonmonotone." (are the fastest and more efficient ones)
			minmover.tolerance(tols[i]) # Set minimization tolerances (0.01 --> 1% from the predicted minimum)
			minmover.apply(pose) # <-- HERE IT MINIMIZES
			score2 = scorefxn(pose)  # Energy after minimization
			tf = time.time()
			print("%d,%d) Energy drops in MinMover from %10.3f to %10.3f in %.3f seconds." % (x, i, score1, score2, tf-t0))
			time_minimization += tf-t0

			#out_name = "min_relax4_%05d.pdb" % (i+1)
			#pose.dump_pdb(out_name)

			#exit()


		this_score = score2
		if this_score < best_score:  # if lower energy...
			best_pose.assign(pose)
			best_score = this_score
			print("%d) Good score: %f" % (x, best_score))
		else:
			pose.assign(best_pose) # if bad, you take the best one...
			if premature and x > 1: # If not first...
				print("%d) Bad score: %f --> Forcing exit (--premature enabled)" % (x, this_score))
				break;
			else:
				print("%d) Bad score: %f" % (x, this_score))

	scorefxn.set_weight(core.scoring.fa_rep,fa_rep_wt)  # restore fa weight
	print("BEST SCORE: %f  (repacking time: %.3f s,  minimization time: %.3f s" % (best_score,time_packing,time_minimization))
	return(best_score)



