r"""
does simulated annealing to find identities for surface residues that 
minimize the fluctuations of average electrostatic potential on the solvent accessible
surface of a symmetric protein multimer. Call as:

.. code-block:: python3

	python chargedesign.py

Inputfile is `input/chargedesign.json`, sample imputfile:

.. code-block:: javascript

	{
	"pdb_in" : "pdb/in/M_C3_Q6.pdb",
	"pdb_out_folder" : "pdb/out/",
	"pdb_out_prefix" : "_",
	"symmfile" : "sym/c3.symm",
	"designable" : 
	[
	257,258,261,264,265,268,269,272,17,21,25,28,29,31,32,34,38,42,62,65,68,
	72,76,79,83,87,90,93,94,96,103,124,125,127,130,134,138,141,145,149,152,
	155,186,187,189,196,200,203,207,211,214,217,218,219,230,237,241,245,248,
	249,250,251,253,254
	],
	"average_potentials_exclude_residues" : 
	[
	39,  40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
	101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,
	163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178
	],
	"chargefile" : "constraints/M_chargeadjust.charge",
	"compfile" : "constraints/M_chargeadjust.comp",
	"scorefile" : "scores/M_chargedesign_scores.txt",
	"kT_sampling" : 0.005,
	"n_mutate" : 4,
	"n_steps" : 10000	
	}

Where:

* `pdb_in` is  a pdb file containing the input structure, which should be a Cn symmetric multimer.
* `pdb_out_folder` is the folder where output structures will be written
* `pdb_out_prefix` is used in the output filename in between the `pdb_in` filename and the output sequence number 
* `symmfile` is the rosetta symmetry file that turns the monomer (chain A of pdb_in) into the correct multimer structure
* `designable` is the set of (surface) residues that can in principle be mutated. This can be generated with the script `getresiduesets.py`
* `average_potentials_exclude_residues` is the residueset to be excluded from the calculation of the averages of the electrostatic potential
* `chargefile` is the rosetta charge file that constrains the net charge of the protein
* `compfile` is the rosetta comp file that constrains the aminoacid composition of the protein
* `scorefile` has the output scores
* `kT_sampling` controls the fluctuations of the score allowed in the simulated annealing search.
* `n_mutate` is the number of mutations per simulated annealing step, residues to be mutated are selected randomly from the `designable` residueset.
* `n_steps` is the number of simulated annealing steps
	 
	 
"""

import pyrosetta
import json
import math
import string
import os
import sys
import numpy
import random

# local
import potentials
import residuesets

# utility
def charge(sequence):
	ch = 0
	for aa in sequence:
		if aa in ['R','K']: ch+=1
		elif aa in ['E','D']: ch-=1
	return ch

def nplus(sequence):
	np = 0
	for aa in sequence:
		if aa in ['R','K']: np+=1
	return np

def nmin(sequence):
	nm = 0
	for aa in sequence:
		if aa in ['D','E']: nm+=1
	return nm

if __name__ == "__main__":

	# main
	scriptname = os.path.basename(__file__).rsplit(".",maxsplit=1)[0]
	thisscript = "mutates to minimize spatial fluctuations of the surface potential, run as python chargedesign.py, input in input/chargedesign.json"
	inputdir = "input"
	inputfile = os.path.join(inputdir,scriptname+".json")

	print("script:", scriptname)
	print(thisscript)
	print("inputfile: ",inputfile)

	print('\n ====> chargedesign setup: read & parse json input file\n')
	with open(inputfile, 'r') as f:
		input = json.load(f)
	pdb_in = input["pdb_in"]
	pdb_out_folder = input["pdb_out_folder"]
	pdb_out_prefix = input["pdb_out_prefix"]
	symmfile = input["symmfile"]
	designable = input["designable"]
	average_potentials_exclude_residues = input["average_potentials_exclude_residues"]
	chargefile = input["chargefile"]
	compfile = input["compfile"]
	scorefile = input["scorefile"] 
	n_steps = input["n_steps"]
	n_mutate = input["n_mutate"]
	kT_sampling = input["kT_sampling"]

	# derived/other vars
	pdb_filename_no_extension = os.path.basename(pdb_in).rsplit(".",maxsplit=1)[0]
	def pdb_out(n):
		return pdb_out_folder+pdb_filename_no_extension+pdb_out_prefix+str(n+1)+".pdb"
	resfile = "res/"+scriptname.rsplit(".",maxsplit=1)[0]+".res"

	# start of script 
	print ("\n ====> chargedesign setup: Initialize pyRosetta \n")
	pyrosetta.init()

	# scorefunctions 
	print("\n ====>  chargedesign setup: Get full Atom Score Function \n")
	scorefxn = pyrosetta.get_fa_scorefxn()
	scorefxn.set_weight(pyrosetta.rosetta.core.scoring.netcharge, 1.0) # for NetChargeConstraintMover
	scorefxn.set_weight(pyrosetta.rosetta.core.scoring.aa_composition, 1.0) # for CompositionContraintMover
	scorefxn_noconstraints = pyrosetta.get_fa_scorefxn()

	print("\n ====>  chargedesign setup: load multimer, extract chain A and symmetrize using symmetry file ")
	pose_in =  pyrosetta.io.pose_from_pdb(pdb_in)
	pose_unit = pyrosetta.rosetta.core.pose.Pose()
	nresidues = pose_unit.total_residue() 
	chain_A_id =  pyrosetta.rosetta.core.pose.get_chain_id_from_chain('A',pose_in)
	pyrosetta.rosetta.core.pose.append_pose_to_pose(pose_unit, pose_in.split_by_chain(chain_A_id))
	pose_symm = pyrosetta.rosetta.core.pose.Pose()
	pose_symm.assign(pose_unit)
	symmetrize = pyrosetta.rosetta.protocols.symmetry.SetupForSymmetryMover(symmfile)
	symmetrize.apply(pose_symm)
	
	# charge constraints - to whole protein
	print("\n ===>  chargedesign setup: Set up NetChargeConstraintMover\n")
	charge_constraint = pyrosetta.rosetta.protocols.aa_composition.AddNetChargeConstraintMover()
	charge_constraint.create_constraint_from_file(chargefile)
	charge_constraint.apply(pose_symm)

	# composition constraints - to whole protein 
	print("\n ===>  chargedesign setup: Set up CompositionConstraintMover\n")
	composition_constraint = pyrosetta.rosetta.protocols.aa_composition.AddCompositionConstraintMover()
	composition_constraint.create_constraint_from_file(compfile)
	composition_constraint.apply(pose_symm)       

	print("\n ===>  chargedesign setup: initialize score variables, save scores of initial structure\n")

	# initialize score variables

	pose_symm.dump_pdb(pdb_out(-1))
	pot = potentials.Potentials(pdb_out(-1),average_potentials_exclude_residues)
	psi_av = pot.psi_av
	delta_psi_sq = pot.delta_psi_sq 
	rosetta_score = scorefxn_noconstraints(pose_symm)
	delta_score = 0.0
	pose_symm_old = pyrosetta.rosetta.core.pose.Pose()

	# save scores for initial structure
	header = "pdb_file		Q		N-		N+		rosetta_score		psi_av		delta_psi_sq\n" 
	f = open(scorefile,'w')
	f.write(header)
	f.write(pdb_out(-1)+"	"+str(charge(pose_symm.sequence()))+"	"+str(nmin(pose_symm.sequence()))+"		"+str(nplus(pose_symm.sequence())))
	f.write("		"+str(rosetta_score)+"		"+str(psi_av)+"		"+str(delta_psi_sq)+"\n")
	f.close()

	print("\n ===>  chargedesign: start simulated annealing steps\n")
	
	# pack rotamers   
	for n in range(n_steps):
	
		# save stuff from previous round
		psi_av_old = psi_av
		delta_psi_sq_old = delta_psi_sq
		pose_symm_old.assign(pose_symm)
		rosetta_score_old = rosetta_score

		print("\n ====> chargedesign annealing step: pick random resiudes to be mutated, find neighbouring residues for repacking, round "+str(n+1)+"\n")
		mutate = list(numpy.random.choice(designable, size=n_mutate, replace=False))
		distance_threshold = 7.0
		close = residuesets.residueset_close_to_target_residueset(pose_unit,mutate,distance_threshold)
		repack = list(set(close).difference(set(mutate)))
		print("mutate = ", mutate)
		print("repack = ", repack)
		residueset_nataa = repack
		residueset_allaa = mutate
		residuesets.residueset_to_resfile(resfile,residueset_nataa,residueset_allaa)
	
		print("\n ====> chargedesign annealing step: setup and apply PackRotamersMover, round "+str(n+1)+"\n")
		task_pack = pyrosetta.standard_packer_task(pose_symm)
		pyrosetta.rosetta.core.pack.task.parse_resfile(pose_symm, task_pack, resfile)
		packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(scorefxn, task_pack)
		packer.apply(pose_symm) 
		pose_symm.dump_pdb(pdb_out(n))
		print("\n...PackRotamersMover done\n")
	
		print("\n ====> chargedesign annealing step: calculate potential averages and rosetta score, round "+str(n+1)+"\n")
		pot = potentials.Potentials(pdb_out(n),average_potentials_exclude_residues)
		psi_av = pot.psi_av
		delta_psi_sq = pot.delta_psi_sq
		rosetta_score = scorefxn_noconstraints(pose_symm)
	
		print("\n ====> chargedesign annealing step: calculate changes in score, round "+str(n+1)+"\n")
		delta_score = delta_psi_sq - delta_psi_sq_old
		boltzmann = math.exp(-delta_score/kT_sampling)
		xi = random.random()
		print("\n===========================")
		print("delta_score = ",delta_score)
		print("kT_sampling = ",kT_sampling)
		print("============================\n")	
	
		print("\n ====> chargedesign annealing step: apply Metropolis acceptance criterium, round "+str(n+1)+"\n") 
		# first step (n=0) always accept to allow for changing the charge 
		if (delta_score < 0.0 or boltzmann > xi or n==0): 
			# accept -save scores to scorefile
			print("\n==================")
			print("step ",n+1,"ACCEPTED")
			print("==================\n")
			with open(scorefile, "a") as f:
				f.write(pdb_out(n)+"	"+str(charge(pose_symm.sequence()))+"	"+str(nmin(pose_symm.sequence()))+"		"+str(nplus(pose_symm.sequence())))
				f.write("		"+str(rosetta_score)+"		"+str(psi_av)+"		"+str(delta_psi_sq)+"\n")
		else: 
			# reject - revert to old values
			pose_symm.assign(pose_symm_old)
			psi_av = psi_av_old
			delta_psi_sq = delta_psi_sq_old
			rosetta_score = rosetta_score_old
			pose_symm.dump_pdb(pdb_out(n))
			print("\n==================")
			print("step ",n+1,"REJECTED")
			print("==================\n")
	
	print("done.") 


