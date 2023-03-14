""" 

module containing set of utility functions used in the main script chargedesign.py

"""

import numpy
import scipy
import gridData
import pyrosetta

import dx
import sas
import potentials
import pdbatom

def check_symmetry_file(pose_in, symm_file):
	"""
	
	checks symmetry file symm_file: should be for Cn symmetry, with n being equal to the number of chains in pose_in 

	:param pose_in: pose to be tested
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
	
	:param symm_file: path to rosetta symmetry file 
	:type symm_file: `str`

	"""	
	with open(symm_file, 'r') as f:
		symm_lines = f.readlines()
	found = False
	for line in symm_lines:
		if line.strip().startswith("symmetry_name") and len(line.split())==2:
			symm = line.split()[1]
			found = True
	if not found: 
		print("...chargedesign_utilities.check_symmetry_file: did not find properly formatted `symmetry_name' line in symm_file.")
		exit(1)
	if not (symm == "c"+str(pose_in.num_chains())):
		print("...chargedesign_utilities.check_symmetry_file: cn symmetry required with n equal to number of chains in pdb_in.")
		exit(1)
	print("...chargedesign_utilities:check_symmetry_file: number of chains in pdb_in  = ",pose_in.num_chains(),", symmetry ", symm)

def get_scores(pose,pdb_out,exclude_residues,scorefxn):
	"""
	
	calculates scores after chargedesign metropolis monte-carlo trial move has been performed. 

	:param pose: pose after the trial move, used to get sequence and rosetta_energy. 
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
	
	:param pdb_out: path to which the pose will be saved before score calculation (needed for pqr2pdb/apbs calculations which start from pdb file) 
	:type pdb_out: `str`

	:param exclude_residues: list of pose residue indices (for asymmetric unit) for residues to be excluded from the calculation of averages of the 	electrostatic potential over the protein solvent acessible surface. 
	:type exclude_residues: `list` of `int`

	:param scorefxn: pyrosetta scorefunction for calculating the rosetta_energy
	:type scorefxn: `pyrosetta.rosetta.core.scoring.ScoreFunction`
	
	"""	

	scores  = dict()
	
	# first save pose to pdb_out, since score computation also requires pdb file
	pose.dump_pdb(pdb_out)
	
	# calculate electrostatic potential averages
	print("...chargedesign_utilities.get_score: calculate averages of electrostatic potential at solvent acessible surface")
	pot = potentials.Potentials(pdb_out,exclude_residues)
	scores["psi_av"] = pot.psi_av
	scores["delta_psi_sq"] = pot.delta_psi_sq 
	
	# calculate rosetta energy
	print("...chargedesign_utilities.get_score: calculate rosetta energy")
	scores["rosetta_energy"] = scorefxn(pose)
	
	# get sequence scores
	print("...chargedesign_utilities.get_score: calculate sequence scores")
	scores["rosetta_energy"] = scorefxn(pose)
	charge = 0
	n_plus = 0
	n_min = 0
	for aa in pose.sequence():
		if aa in ['R','K']: 
			charge+=1
			n_plus+=1
		elif aa in ['E','D']: 
			charge-=1
			n_min+=1
	scores["charge"] = charge
	scores["n_plus"] = n_plus
	scores["n_min"] = n_min
	
	return scores
	
def save_scores(score_file,scores,pdb_out,n):
	"""
	
	saves scores of chargedesign metropolis monte-carlo run to score_file. 

	:param score_file: path to score_file 
	:type score_file: `str`
	
	:param scores: dictionary with score information to be saved  
	:type scores: `dict`

	:param pdb_out: path to pdb_file to which output of the metropolis monte-carlo trial move (step n) was saved.  
	:type pdb_out: `str`

	:param n: step number of the metropolis monte-carlo run 
	:type n: `int`
	
	"""	
	
	# at start create file and write header, otherwise append
	if n == 0:
		header = "step		pdb_file		charge		n_min		n_plus		rosetta_energy (kT)		psi_av (mV)		delta_psi_sq (mV^2)\n" 
		f = open(score_file,'w')
		f.write(header)
	else:
		f = open(score_file,'a')
	
	# write scores and close file (to avoid losing data)
	f.write(str(n)+"   ")
	f.write(pdb_out+"   ")
	f.write(str(scores["charge"])+"   ")
	f.write(str(scores["n_min"])+"   ")
	f.write(str(scores["n_plus"])+"   ")
	f.write(str(scores["rosetta_energy"])+"   ")
	f.write(str(scores["psi_av"])+"   ")
	f.write(str(scores["delta_psi_sq"])+"\n")
	f.close()

def find_neighbours(pose,target_residueset,distance_threshold):
	"""
	
	find residueset of residues close to residues in a target_residueset given some distance_treshold. 
	The target_residueset itself is excluded from the neighbour list. 

	:param pose: input pose 
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
	
	:param target_residueset: residueset of residues close to which we are seeking other residues 
	:type target_residueset: `list` of `int`

	:param distance_treshold: distance_treshold in Angstrom 
	:type distance_treshold: `float`

	:return: residueset
	:rtype: `list` of `int`
	
	"""	
	# make index selector for target_residueset
	n_res = pose.total_residue() 
	index_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
	index_string = ",".join([str(n) for n in target_residueset])
	index_selector.set_index(index_string)
	
	# make contact selector, use index_selector to set target_residueset as central group
	contact_selector = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
	contact_selector.threshold(distance_threshold)
	contact_selector.central_residue_group_selector(index_selector)
	
	# create list of contact residues, exclude target_residueset from the list.
	contact_residues = contact_selector.apply(pose)
	contact_residueset = set()
	for n in range(n_res):
		if contact_residues[n+1]: contact_residueset.add(n+1)
	return list(contact_residueset.difference(target_residueset))

def write_resfile(filename,residueset_nataa, residueset_allaa):
	"""
	creates rosetta resfile, given residuesets with NATAA and ALAA entries
	
	:param filename: path for resfile 
	:type filename: `str`
	
	:param residueset_nataa: residueset with NATAA entries 
	:type residueset_nataa: `list` of `int`
	
	:param residueset_allaa: residueset with ALLAA entries 
	:type residueset_allaa: `list` of `int`
	
	"""
	f = open(filename,'w')
	f.write("NATRO\n")
	f.write("USE_INPUT_SC\n")
	f.write("start\n")
	for index in residueset_nataa:
		f.write(str(index)+"   A    NATAA\n")
	for index in residueset_allaa:
		f.write(str(index)+"   A    ALLAA\n")
	f.close()