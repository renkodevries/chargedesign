""" 

module containing class ChargeDesigner

"""

#standard
import json
import os
import math
import random

# external
import numpy
import scipy
import gridData
import pyrosetta

# local 
import dx
import sas
import potentials
import pdbatom


class ChargeDesigner:
	r"""
		
		Class for doing charge design: minimization of fluctuations of electrostatic potential on protein solvent acessible surface
		
		:Usage example:
	
		.. code-block:: python3
		
			# basic usage
			chargedesigner = chargedesigner.ChargeDesigner(input_file)
			while chargedesigner.step < chargedesigner.n_total_steps:
				chargedesigner.step+=1
				chargedesigner.set_kT_sampling()
				chargedesigner.do_trial_move()
				chargedesigner.set_scores()
				chargedesigner.apply_acceptance_criterium()
				chargedesigner.save_scores()
				chargedesigner.check_for_new_minimum()

	
	
		:Sample input file: 
		
		.. code-block:: javascript

			{
			"pdb_in" : "pdb/in/M_C3.pdb",
			"pdb_out_folder" : "pdb/out/",
			"run_id" : "M_C3_Q12",
			"symm_file" : "sym/c3.symm",
			"charge_file" : "constraints/M_chargeadjust.charge",
			"comp_file" : "constraints/M_chargeadjust.comp",
			"score_file" : "scores/M_chargedesign_scores.txt",
			"designable_residues" : 
			[
			257,258,261,264,265,268,269,272,17,21,25,28,29,31,32,34,38,42,62,65,68,
			72,76,79,83,87,90,93,94,96,103,124,125,127,130,134,138,141,145,149,152,
			155,186,187,189,196,200,203,207,211,214,217,218,219,230,237,241,245,248,
			249,250,251,253,254
			],
			"exclude_residues" : 
			[
			39,  40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54,
			101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,
			163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178
			],
			"repack_distance_treshold" : 7.0,
			"n_mutate" : 4,
			"kT_sampling_start" : 1,
			"kT_sampling_end" : 0.1,
			"n_steps_hold_start" : 100,
			"n_steps_hold_end" : 1000,
			"n_steps_linear_gradient" : 100
			}

		Where:

		* `pdb_in` is  a pdb file containing the input structure, which should be a Cn symmetric multimer.
		* `pdb_out_folder` is the folder where output structure will be written, also used for writing and reading temporary pdb files
		* `run_id` is used for naming various files, should be unique if multiple copies of the program run side-by-side. Filenames: output filename is run_id.pdb, temporary res_file and pdb_file step 5 are run_id_5_temp.res and run_id_5_temp.pdb 
		* `symm_file` is the rosetta symmetry file that turns the monomer (first chain of pdb_in) into the correct multimer structure
		* `charge_file` is the rosetta charge file that constrains the net charge of the protein
		* `comp_file` is the rosetta comp file that constrains the aminoacid composition of the protein
		* `score_file` has the output scores
		* `designable_residues` is the set of (surface) residues that can in principle be mutated.
		* `exclude_residues` is the residueset to be excluded from the calculation of the averages of the electrostatic potential
		* `repack_distance_treshold` residues within this distance (in Angstrom) from the residues to be mutated, will be repacked.
		* `n_mutate` is the number of mutations per simulated annealing step, residues to be mutated are selected randomly from the `designable` residueset.
		* `kT_sampling_start` . kT_sampling controls the fluctuations of the score allowed in the simulated annealing search, this is the starting value.
		* `kT_sampling_end` final value of kT_sampling.
		* `n_steps_hold_start`  number of initial steps at `kT_sampling_start`
		* `n_steps_hold_end`  number of final steps at `kT_sampling_end`
		* `n_steps_linear_gradient`  number of steps of linear gradient from `kT_sampling_start` to `kT_sampling_end`
	
	"""
	def __init__(self,input_file):
		"""
		Constructor, initializes from json input file
		
		:param input_file: name of json input file  
		:type input_file: `str`
	
		"""
		print("\n...chargedesigner.__init__: reading input, initializing variables\n")
		with open(input_file, 'r') as f:
			self.input = json.load(f)
		self.pdb_in = self.input["pdb_in"]
		self.pdb_out_folder = self.input["pdb_out_folder"]
		self.run_id = self.input["run_id"]
		self.symm_file = self.input["symm_file"]
		self.charge_file = self.input["charge_file"]
		self.comp_file = self.input["comp_file"]
		self.score_file = self.input["score_file"] 
		self.designable_residues = self.input["designable_residues"]
		self.exclude_residues = self.input["exclude_residues"]
		self.repack_distance_threshold = self.input["repack_distance_treshold"]
		self.n_mutate = self.input["n_mutate"]
		self.n_steps_hold_start = self.input["n_steps_hold_start"]
		self.n_steps_hold_end = self.input["n_steps_hold_end"]
		self.n_steps_linear_gradient = self.input["n_steps_linear_gradient"]
		self.kT_sampling_start = self.input["kT_sampling_start"]
		self.kT_sampling_end = self.input["kT_sampling_end"]
		self.n_total_steps = self.n_steps_hold_start + self.n_steps_linear_gradient + self.n_steps_hold_end
		self.step = 0
		self.kT_sampling = self.kT_sampling_start
		self.mutate = []
		self.repack = []
		self.pose = pyrosetta.rosetta.core.pose.Pose()
		self.pose_old = pyrosetta.rosetta.core.pose.Pose()
		self.scores = dict()
		self.scores_old = dict()
		self.scores_min = dict()
		self.init_scores()
		self.scores_old = self.scores.copy()
		self.scores_min = self.scores.copy()
		self.n_accepted = 0		
		print("\n...chargedesigner.__init__: getting rosetta scorefunctions\n")
		self.scorefxn = pyrosetta.get_fa_scorefxn()
		self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.netcharge, 1.0) # for NetChargeConstraintMover
		self.scorefxn.set_weight(pyrosetta.rosetta.core.scoring.aa_composition, 1.0) # for CompositionContraintMover
		self.scorefxn_noconstraints = pyrosetta.get_fa_scorefxn()
		print("\n...chargedesigner.__init__: getting asymmetric unit from pdb_in\n")
		self.pose_in = pyrosetta.io.pose_from_pdb(self.pdb_in)
		self.n_chains = self.pose_in.num_chains()
		self.pose_unit = pyrosetta.rosetta.core.pose.Pose()
		pyrosetta.rosetta.core.pose.append_pose_to_pose(self.pose_unit, self.pose_in.split_by_chain(1))
		self.n_res = self.pose_unit.total_residue()
		print("\n...chargedesigner.__init__: symmetrize asymmetric unit using symmetry file symm_in\n")
		self.pose.assign(self.pose_unit)
		self.check_symmetry_file()	
		self.symmetrize_pose()
		print("\n...chargedesigner.__init__: apply charge and composition constraints to scoring for symmetric pose\n")
		self.apply_charge_constraint_to_pose()
		self.apply_composition_constraint_to_pose()
	# ===============
	# HELPER METHODS
	# ===============
	def init_scores(self):
		self.scores["step"] = 0
		self.scores["charge"] = 0
		self.scores["n_min"] = 0
		self.scores["n_plus"] = 0
		self.scores["rosetta_energy"] = 0
		self.scores["psi_av"] = 0
		self.scores["delta_psi_sq"] = 0
		self.scores["kT_sampling"] = 0
		self.scores["n_accepted"] = 0
	def apply_charge_constraint_to_pose(self):
		charge_constraint = pyrosetta.rosetta.protocols.aa_composition.AddNetChargeConstraintMover()
		charge_constraint.create_constraint_from_file(self.charge_file)
		charge_constraint.apply(self.pose)
	def apply_composition_constraint_to_pose(self):
		composition_constraint = pyrosetta.rosetta.protocols.aa_composition.AddCompositionConstraintMover()
		composition_constraint.create_constraint_from_file(self.comp_file)
		composition_constraint.apply(self.pose)		  
	def symmetrize_pose(self):
		symmetrize = pyrosetta.rosetta.protocols.symmetry.SetupForSymmetryMover(self.symm_file)
		symmetrize.apply(self.pose)
	def check_symmetry_file(self):
		with open(self.symm_file, 'r') as f:
			symm_lines = f.readlines()
		found = False
		for line in symm_lines:
			if line.strip().startswith("symmetry_name") and len(line.split())==2:
				symm = line.split()[1]
				found = True
		if not found: 
			print("...chargedesigner.check_symmetry_file: did not find properly formatted `symmetry_name' line in symm_file.")
			exit(1)
		if not (symm == "c"+str(self.n_chains)):
			print("...chargedesigner.check_symmetry_file: cn symmetry required with n equal to number of chains in pdb_in.")
			exit(1)
		print("\n...chargedesigner.check_symmetry_file: number of chains in pdb_in	= ",self.n_chains,", symmetry ", symm,"\n")
	def find_neighbours(self,target_residueset):
		# make index selector for target_residueset
		index_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
		index_string = ",".join([str(n) for n in target_residueset])
		index_selector.set_index(index_string)
		# make contact selector, use index_selector to set target_residueset as central group
		contact_selector = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
		contact_selector.threshold(self.repack_distance_threshold)
		contact_selector.central_residue_group_selector(index_selector)
		# create list of contact residues, exclude target_residueset from the list.
		contact_residues = contact_selector.apply(self.pose_unit)
		contact_residueset = set()
		for n in range(self.n_res):
			if contact_residues[n+1]: contact_residueset.add(n+1)
		return list(contact_residueset.difference(target_residueset))
	def write_resfile(self, res_file, residueset_nataa, residueset_allaa):
		f = open(res_file,'w')
		f.write("NATRO\n")
		f.write("USE_INPUT_SC\n")
		f.write("start\n")
		for index in residueset_nataa:
			f.write(str(index)+"   A	NATAA\n")
		for index in residueset_allaa:
			f.write(str(index)+"   A	ALLAA\n")
		f.close()
	# ===============
	#  MAIN METHODS
	# ===============
	def set_kT_sampling(self):
		"""
		
		Sets the actual kT_sampling, depending on value of step
		
		"""
		delta_kT = self.kT_sampling_start - self.kT_sampling_end
		kT_step = delta_kT/float(self.n_steps_linear_gradient)
		dn = float(self.step - self.n_steps_hold_start)
		if self.step <= self.n_steps_hold_start:
			self.kT_sampling = self.kT_sampling_start
		elif self.step <= self.n_steps_hold_start + self.n_steps_linear_gradient:
			self.kT_sampling = self.kT_sampling_start - dn*kT_step
		else:
			self.kT_sampling = self.kT_sampling_end
		print("\n...chargedesigner.set_kT_sampling: kT_sampling = ",self.kT_sampling)
	def do_trial_move(self):
		"""
		
		Does chargedesign trial move: picks n_mutate residues from the residueset designable_residues (residueset "mutate"),
		finds the residues within a distance repack_distance_threshold from the residues being mutated (residueset "repack"),
		writes a temporary resfile (based on the "mutate" and "repack" residuesets), sets up a PackRotamersMover, 
		and applies this to the pose. The pose is saved to a temporary pdb file, the resfile is deleted.
		
		"""
		# temp file names	
		res_dir = "res"
		res_file = os.path.join(res_dir,self.run_id+"_"+str(self.step)+"_temp.res")
		pdb_file = os.path.join(self.pdb_out_folder,self.run_id+"_"+str(self.step)+"_temp.pdb")
		# action
		mutate = list(numpy.random.choice(self.designable_residues, size=self.n_mutate, replace=False))
		repack = self.find_neighbours(mutate)
		print("\n...mutate = ", mutate)
		print("...repack = ", repack,"\n")
		self.write_resfile(res_file,repack,mutate)
		print("...saved temporary res file ",res_file,"\n")
		task_pack = pyrosetta.standard_packer_task(self.pose)
		pyrosetta.rosetta.core.pack.task.parse_resfile(self.pose, task_pack, res_file)
		packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover(self.scorefxn, task_pack)
		packer.apply(self.pose) 
		self.pose.dump_pdb(pdb_file)
		os.remove(res_file)
		print("\n...saved temporary pdb file ",pdb_file)
		print("...removed temporary res file ",res_file)
		print("...PackRotamersMover done\n")
	def set_scores(self):
		"""
		
		Calculates the scores to be optimized as well as any other relevant scores. 
		Uses class potentials.Potentials to get the averages of the electrostatic 
		potential on the solvent accessible surface of the protein. After scores are calculated,
		temporary pdb file for pose (used by potentials.Potentials) is removed.
		
		"""
		# temp pdb file for calculating electrostatic potential
		pdb_file = os.path.join(self.pdb_out_folder,self.run_id+"_"+str(self.step)+"_temp.pdb")
		# calculate electrostatic potential averages
		print("\n...chargedesigner.set_scores: calculate scores from electrostatic potential\n")
		pot = potentials.Potentials(pdb_file,self.exclude_residues)
		self.scores["psi_av"] = pot.psi_av
		self.scores["delta_psi_sq"] = pot.delta_psi_sq 
		# calculate rosetta energy
		print("...chargedesigner.set_scores: calculate rosetta energy")
		self.scores["rosetta_energy"] = self.scorefxn(self.pose)
		# get sequence scores
		print("...chargedesigner.set_scores: calculate sequence scores")
		charge = 0
		n_plus = 0
		n_min = 0
		for aa in self.pose.sequence():
			if aa in ['R','K']: 
				charge+=1
				n_plus+=1
			elif aa in ['E','D']: 
				charge-=1
				n_min+=1
		self.scores["charge"] = charge
		self.scores["n_plus"] = n_plus
		self.scores["n_min"] = n_min
		self.scores["kT_sampling"] = self.kT_sampling
		self.scores["step"] = self.step
		self.scores["n_accepted"] = self.n_accepted
		# temp pdb file no longer needed
		os.remove(pdb_file)
		print("...removed temporary pdb file ",pdb_file)
	def apply_acceptance_criterium(self):
		"""
		
		Applies Metropolis acceptance criterium to the new pose (after do_trial_move and set_scores have been called).
		The target score to be minimized is delta_psi_sq, the mean square of the fluctuations of the elecrostatic potential 
		on the solvent accessible surface of the protein. If rejected, the trial move is made undone.
		
		"""
		delta = self.scores["delta_psi_sq"] - self.scores_old["delta_psi_sq"]
		boltzmann = math.exp(-delta/self.kT_sampling)
		xi = random.random()
		if (self.step==1):
			print("...chargedesigner.apply_acceptance_criterium: FIRST STEP, ALWYAS ACCEPTED\n")
		elif (delta < 0.0 or boltzmann > xi): 
			print("...chargedesigner.apply_acceptance_criterium, delta_psi_sq (old) = ",self.scores_old["delta_psi_sq"])
			print("...chargedesigner.apply_acceptance_criterium, delta_psi_sq (new) = ",self.scores["delta_psi_sq"])
			print("...chargedesigner.apply_acceptance_criterium, delta/kT_sampling = ",delta/self.kT_sampling)
			print("...chargedesigner.apply_acceptance_criterium: TRIAL MOVE ACCEPTED\n")
			self.n_accepted+=1
		else: 
			print("...chargedesigner.apply_acceptance_criterium, delta_psi_sq (old) = ",self.scores_old["delta_psi_sq"])
			print("...chargedesigner.apply_acceptance_criterium, delta_psi_sq (new) = ",self.scores["delta_psi_sq"])
			print("...chargedesigner.apply_acceptance_criterium, delta/kT_sampling = ",delta/self.kT_sampling)
			print("...chargedesigner.apply_acceptance_criterium: TRIAL MOVE REJECTED\n")
			self.pose.assign(self.pose_old)
			self.scores = self.scores_old.copy()
	def save_scores(self):
		"""
		
		Saves scores to score_file. Saves current pose and scores to be used
		as "old" values for the next step.
		
		"""
		# at start create file and write header, otherwise append
		step_str = '{:12d}'.format(self.step) 
		charge_str = '{:12d}'.format(self.scores["charge"]) 
		n_plus_str = '{:12d}'.format(self.scores["n_plus"]) 
		n_min_str = '{:12d}'.format(self.scores["n_min"]) 
		rosetta_energy_str = '{:24.2f}'.format(self.scores["rosetta_energy"]) 
		psi_av_str = '{:16.2f}'.format(self.scores["psi_av"])
		delta_psi_sq_str = '{:24.2f}'.format(self.scores["delta_psi_sq"])
		kT_sampling_str = '{:16.2f}'.format(self.scores["kT_sampling"])
		n_accepted_str = '{:16d}'.format(self.scores["n_accepted"])
		dataline = step_str+charge_str+n_plus_str+n_min_str+rosetta_energy_str+psi_av_str+delta_psi_sq_str+kT_sampling_str+n_accepted_str
		if self.step == 1:
			header = "step".rjust(12)
			header += "charge".rjust(12)
			header += "n_plus".rjust(12)
			header += "n_min".rjust(12)
			header += "rosetta_energy(kT)".rjust(24)
			header += "psi_av(mV)".rjust(16)
			header += "delta_psi_sq(mV^2)".rjust(24)
			header += "kT_sampling".rjust(16)
			header += "n_accepted".rjust(16)
			f = open(self.score_file,'w')
			f.write(header+"\n")
			f.write(dataline+"\n")
		else:
			f = open(self.score_file,'a')
			f.write(dataline+"\n")
			f.close()
		# prepare for next step
		self.scores_old = self.scores.copy()
		self.pose_old.assign(self.pose)
	
	def check_for_new_minimum(self):
		"""
		
		Checks if the target score (delta_psi_sq) has a new lowest value. If so, 
		the new minimum is saved and the new minimum structure is saved to disk as pdb file.
		
		"""
		pdb_out = os.path.join(self.pdb_out_folder,self.run_id+".pdb")
		if (self.step==1):
			# delete if pdb with same new still exists from previous run
			if os.path.isfile(pdb_out): os.remove(pdb_out)
			#  step==1 pose is first minimum structure
			self.pose.dump_pdb(pdb_out)
			self.scores_min = self.scores.copy()
		elif self.scores["delta_psi_sq"] <  self.scores_min["delta_psi_sq"]:
			# new min: remove old structure, replace by new min
			os.remove(pdb_out)
			self.pose.dump_pdb(pdb_out)
			self.scores_min = self.scores.copy()
			print("...chargedesigner.check_for_new_minimum: new minimum, delta_psi_sq = ",self.scores["delta_psi_sq"],"mV^2")
		
		
		