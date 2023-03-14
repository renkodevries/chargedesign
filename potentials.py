""" 

module containing class Potentials

"""

import gridData 
import numpy
import scipy.interpolate
import math
import subprocess
import os
import pdbatom
import dx
import sas

class Potentials:
	"""Class for calculating averages of the electrostatic potentials of proteins at 
	their solvent accessible surfaces. Computes the average potential `psi_av` and the
	its mean square fluctuations `delta_psi_sq`. Calls pqr2pdb to convert the input pdb file 
	into pqr file with atom charge. Next calls APBS which takes the pqr file as input and
	produces a .dx file with values of the electrostatic potential on a regular spatial grid.
	Next calculates the averages, using modules dx and sas. Finally deletes all intermediate files 
	no longer needed.

	Requires that pdb2pqr and apbs have been installed. Other (non standard) Python dependencies: 
	gridData, scipy.interpolate, numpy, plus the local modules pdbatom, dx and sas.
	
	:param pdb_in: path to pdb file with input structure
	:type pdb_in: `str`
	
	:param pqr_file: path to pqr file produced by pdb2pqr
	:type pqr_file: `str`
	
	:param dx_file: path to dx file produced by APBS
	:type dx_file: `str`
	
	:param apbs_inputfile: path to APBS input file produced by filling the `apbs_template`
	:type apbs_inputfile: `str`
	
	:param exclude_residues: averages `psi_av` and `delta_psi_sq` are obatined by averaging over all residues minus the exclude_residues 	
	:type exclude_residues: `list` of `int`
	
	:param apbs_template: path to template input file for APBS, this file must already be present.
	:type apbs_template: `str`
	
	:param sas: class instance of sas.Sas for computing solvent accessible surface
	:type sas: `sas.Sas`
	
	:param dx: class instance of dx.Dx for interpolating potentials from APBS.
	:type dx: `dx.Dx`
	
	:param areas: solvent accessible surface area per atom (only heavy atoms)
	:type areas: `list` (length `sas.natom`) of `float`
	
	:param potentials: average potentials at solvent accessible surface per atom (only heavy atoms)
	:type potentials: `list` (length `sas.natom`) of `float`
	
	:param psi_av: average potential,averaged over all residues minus those in exclude_residues
	:type psi_av: `float`
	
	:param delta_psi_sq: mean square potential,averaged over all residues minus those in exclude_residues
	:type delta_psi_sq: `float`
	
	:Example:
	
	.. code-block:: python3
			
		pot = potentials.Potentials("mystruct.pdb",exclude_residues)
		print("psi_av  = ",pot.psi_av)
		print("delta_psi_sq  = ",pot.delta_psi_sq)
				
	"""	
	
	def __init__(self,pdb_file,exclude_residues):
		"""
		Constructor. All computations are done immediately after calling the constructor.
		
		:param pdb_file: path to pdb file with input structure
		:type pdb_file: `str`
		
		:param exclude_residues: list of residue indices (residue sequence numbers from pdb file) 
		for residues to be excluded from computing potential averages 
		:type exclude_residues: `list` of `int`
		
		"""
		self.pdb_file = pdb_file
		self.pdb_file_basename = pdb_file.rsplit(".",maxsplit=1)[0]
		if (len(self.pdb_file_basename.rsplit("/"))>1):
			self.pdb_file_basename = self.pdb_file_basename.rsplit("/")[-1]
		self.pqr_file = "apbs/"+self.pdb_file_basename+".pqr"
		self.dx_file =  "apbs/"+self.pdb_file_basename
		self.apbs_inputfile = "apbs/"+self.pdb_file_basename+"_apbs.in"
		self.exclude_residues = exclude_residues
		self.apbs_template = "apbs/apbs_template.in"
		self.sas = None
		self.dx = None
		self.areas = []
		self.potentials = []
		self.psi_av = float()
		self.delta_psi_sq = float()
		# do calculations
		print("...chargedesign.potentials: clean pdb")
		self.clean_pdb()
		print("...chargedesign.potentials: run_pqr2pdb")
		self.run_pqr2pdb()
		print("...chargedesign.potentials: make_apbs_inputfile")
		self.make_apbs_inputfile()
		print("...chargedesign.potentials: run_apbs")
		self.run_apbs()
		print("...chargedesign.potentials: open apbs map file (dx)")
		self.dx = dx.Dx(self.dx_file+".dx")
		print("...chargedesign.potentials: calculate solvent accessible surface")
		self.sas = sas.Sas(self.pdb_file)
		print("...chargedesign.potentials: calculate average atom potentials")
		self.get_per_atom_areas_and_potentials()
		print("...chargedesign.potentials: calculate average and fluctuations of protein potential")
		self.get_average_and_fluctuations_of_potential()
		print("...chargedesign.potentials: done, removing intermediate files\n")   
		self.remove_intermediate_files()
	def clean_pdb(self):
		"""
		removes any non-pdb lines from the pdb file before running pqr2pdb. 
		Such lines are added by pyrosetta and give warnings in pqr2pdb
		"""
		pdb_keywords  = ["ATOM","EXPDTA","HEADER","TER","REMARK"]
		g = open(self.pdb_file,'r')
		pdb_lines = g.readlines()
		g.close()
		pdb_lines_cleaned=[]
		for line in pdb_lines:
			line_ok = False
			for key in pdb_keywords:
				if line.startswith(key): line_ok = True
			if line_ok: pdb_lines_cleaned.append(line)
		g = open(self.pdb_file,'w')
		for line in pdb_lines_cleaned:
			g.write(line)
		g.close()	
	def get_per_atom_areas_and_potentials(self):
		"""
		computes and sets per atom averages `areas` and `potentials`
		"""
		self.areas = []
		self.potentials = []
		for n in range(self.sas.natoms):
			# area
			selected_dot_indices = self.sas.indices_selected_dots_per_atom[n]
			n_selected_dots =  len(selected_dot_indices)
			dot_positions = self.sas.atom_dots[n]
			whole_atom_sasa = 4*math.pi*(self.sas.atom_radius[n]+self.sas.probe_radius)**2
			area_per_dot = whole_atom_sasa/float(self.sas.dots_per_atom)
			atom_sasa = n_selected_dots*area_per_dot
			self.areas.append(atom_sasa)
			# potential
			selected_dots = self.sas.atom_dots[n][selected_dot_indices,:]
			if n_selected_dots > 0:
				pot = self.dx.value(selected_dots)
				self.potentials.append(pot.sum()/float(n_selected_dots))
			else:
				self.potentials.append(0.0)		
	def get_average_and_fluctuations_of_potential(self):
		"""
		computes and sets overall averages `psi_av` and `delta_psi_sq`
		"""
		# average potential
		psi_av = 0.0
		sum_psi_times_area = 0.0
		sum_area = 0.0
		for n in range(self.sas.natoms):
			at = self.sas.atoms[n]
			ires = at.residueSequenceNumber
			if not ires in self.exclude_residues:
				sum_psi_times_area += self.areas[n]*self.potentials[n]
				sum_area += self.areas[n]
		psi_av = sum_psi_times_area/sum_area
		# fluctuations of potential
		delta_psi_sq = 0.0
		sum_delta_psi_sq_times_area = 0.0
		sum_area = 0.0
		for n in range(self.sas.natoms):
			at = self.sas.atoms[n]
			ires = at.residueSequenceNumber
			if not ires in self.exclude_residues:
				sum_delta_psi_sq_times_area += self.areas[n]*(self.potentials[n]-psi_av)**2
				sum_area += self.areas[n]
		delta_psi_sq = sum_delta_psi_sq_times_area/sum_area
		self.psi_av = psi_av
		self.delta_psi_sq = delta_psi_sq
	def run_pqr2pdb(self):
		"""
		runs pqr2pdb, using same options as used in PyMol `
		"""
		option1 = "--ff=AMBER"
		option2 = "--titration-state-method=propka"
		option3 = "--with-ph=7"
		cmd = "pdb2pqr30"
		do_pqr = subprocess.run([cmd,option1,option2,option3,self.pdb_file,self.pqr_file],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		exitcode = do_pqr.returncode
		if not (exitcode==0):
			print("potentials.run_pqr2pdb: error running pdb2pqr, exitcode "+str(exitcode))
			exit(1)
	def run_apbs(self):
		"""
		runs APBS, using input file `apbs_inputfile`, produced from `apbs_template` `
		"""
		cmd = "apbs"
		do_apbs = subprocess.run([cmd,self.apbs_inputfile],stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
		exitcode = do_apbs.returncode
		if not (exitcode==0):
			print("potentials.run_apbs: error running apbs, exitcode "+str(exitcode))
			exit(1)
	def make_apbs_inputfile(self):
		"""
		produces `apbs_inputfile` from `apbs_template` plus the data items `pqr_file` and `dx_file`
		"""
		# data
		data = dict()
		data['PQR_FILE']= self.pqr_file
		data['DX_FILE']= self.dx_file
		# get template file
		f = open(self.apbs_template,'r')
		template_lines = f.readlines()
		f.close()
		# fill template with data
		for n in range(len(template_lines)):
			line = template_lines[n]
			new_line = line
			for key in data:
				if not line.find(key)==-1:
					replacement = str(data[key])
					new_line = new_line.replace(key,replacement)
			template_lines[n] = new_line 
		# write input file
		f = open(self.apbs_inputfile,'w')
		for line in template_lines:
			f.write(line)
		f.close()
	def remove_intermediate_files(self):
		"""
		removes intermediate files `dx_file` , `pqr_file` , `àpbs_in` , `log_file` 
		"""	
		pqr_f = self.pqr_file
		dx_f = self.dx_file+".dx"
		apbs_in_f = self.apbs_inputfile
		log_f = "apbs/"+self.pdb_file_basename+".log"
		intermediate_files = [pqr_f,dx_f,apbs_in_f,log_f]
		for intermediate_file in intermediate_files:
			if os.path.isfile(intermediate_file): 
				os.remove(intermediate_file)
				print(".......removed ",intermediate_file)
		