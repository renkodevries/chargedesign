""" 

module containing the class Sas

"""

import pdbatom
import numpy
import math
import scipy 

class Sas:
	"""Class for computing solvent accessible surface from pdb files. Dependencies: pdbatom, numpy, math, scipy.
	The solvent accessible surface is calculated as a list of dots. First, for each atom, an exclusion sphere around the atom (with radius equal to the
	vanderwaalsradius of the atom plus the probe radius), is covered with dots. Next, dots of the atom which are occluded by
	exclusion spheres of other atoms are eliminated. 
	
	:param pdb_in: path to pdb file with input structure
	:type pdb_in: `str`
	
	:param atoms: list of atoms as read from pdb_in
	:type atoms: `list` of `PdbAtom`, length `n_atoms`.
	
	:param n_atoms: number of atoms read from pdb_in
	:type n_atoms: `int`
	
	:param atom_radius: list of atom vanderwaals radii
	:type atom_radius: `list` of `float`, length `n_atoms`.
	
	:param atom_position: array holding atom coordinates
	:type atom_position: `numpy.array` with shape (`n_atoms`,3)
	
	:param element_radius: dictionary with atom vanderwaalsradii (values, `float`) of chemical elements (keys, `str`)
	:type element_radius: `dict()`
	
	:param neighbourlist_cutoff: for finding neighbouring atoms (Angstrom)
	:type neighbourlist_cutoff: `float`
	
	:param probe_radius: radius of rolling sphere defining solvent accessible surface area (Angstrom)
	:type probe_radius: `float`
	
	:param dots_per_atom: number of dots used to approximate solvent acessible surface around each atom
	:type dots_per_atom: `int`
	
	:param neighbourlist: for each atom a list of atom indices of neighbouring atoms (within neighbourlist_cutoff)
	:type neighbourlist: `list` of `list` of `int`
	
	:param atom_dots: coordinates of all dots, per atom
	:type atom_dots: `list` of `numpy.array` (list length `n_atoms`), each with shape (`dots_per_atom`,3)
	
	:param indices_selected_dots_per_atom: indices for dots (from atom_dots) which are solvent acessible, per atom
	:type indices_selected_dots_per_atom: `list` (length `n_atoms`) of `list` of `int`
	
	:param self.ndots_sas: total number of dots forming solvent acessible surface
	:type self.ndots_sas: `int`
	
	:Example:
	
	.. code-block:: python3
		
		# lets get the dots forming the surface accesible area of atom 3	
		sas = sas.Sas("trimer.pdb")
		all_dots = sas.atom_dots[3]
		selected_dot_indices = sas.indices_selected_dots_per_atom[3]
		n_surface_dots = len(selected_dot_indices)
		surface_dots = [all_dots[i] for i in selected_dot_indices] # list of numpy arrays with shape (3,) 	
		
	"""	
	def __init__(self,pdb_in):
		"""
		Constructor. All computations are done immediately after calling the constructor.
		
		:param pdb_in: path to pdb file with input structure
		:type pdb_in: `str`

		"""
		# data
		self.pdb_in = pdb_in
		self.atoms = []
		self.natoms = int()
		self.atom_radius = None
		self.atom_position = None
		self.element_radius = dict()
		self.neighbourlist_cutoff = float()
		self.probe_radius = float()
		self.dots_per_atom = int()
		self.neighbourlist = []
		self.atom_dots = []
		self.indices_selected_dots_per_atom = []
		self.ndots_sas = int()
		# initialize
		self.set_default_parameter_values()
		self.load_atoms()
		self.set_atom_position_and_radius_arrays()
		self.build_neighbour_list()
		self.compute_sas()
	def set_default_parameter_values(self):
		r"""
		Sets these default parameter values: 
		
		* `probe_radius` = 1.4 Angstrom 
		* `dots_per_atom` = 50
		* `neighbourlist_cutoff` = 7.0 Angstrom 
		* `element_radius` ['C'] = 1.9 Angstrom
		* `element_radius` ['N']= 1.7 Angstrom    
		* `element_radius` ['O'] = 1.5 Angstrom
		* `element_radius` ['S'] = 1.9 Angstrom
		
		"""
		# lengths in Angstrom
		self.probe_radius = 1.4
		self.dots_per_atom = 50
		self.neighbourlist_cutoff = 7.0
		self.element_radius['C'] = 1.9
		self.element_radius['N'] = 1.7
		self.element_radius['O'] = 1.5
		self.element_radius['S'] = 1.9
	def print_parameter_values(self):
		"""
		Print values of `probe_radius` , `dots_per_atom`,
		`neighbourlist_cutoff` and the `element_radius` dict.
		"""
		print("Sas parameter values:")
		print("probe_radius: ",self.probe_radius)
		print("dots_per_atom: ",self.dots_per_atom)
		print("neighbourlist_cutoff: ", self.neighbourlist_cutoff)
		print("element_radius: ",self.element_radius)
	def set_atom_position_and_radius_arrays(self):
		"""
		Sets `atom_position` and `atom_radius` arrays, called 
		after `load_atoms`.
		"""
		self.atom_position = numpy.zeros((self.natoms,3))
		self.atom_radius = numpy.zeros(self.natoms)
		for n in range(self.natoms):
			at = self.atoms[n]
			x = at.xOrthogonalCoordinate
			y = at.yOrthogonalCoordinate
			z = at.zOrthogonalCoordinate
			el = at.elementSymbol
			self.atom_position[n,:] =numpy.array([x,y,z])
			self.atom_radius[n] = self.element_radius[el]
	def build_neighbour_list(self):		
		"""
		Builds neighbour list for each atom. For each central atom,
		only atoms in its neighbourlist are use to identify dots that
		are solvent accessible.
		"""
		kdtree = scipy.spatial.KDTree(self.atom_position)
		pairs = kdtree.query_pairs(self.neighbourlist_cutoff)
		self.neighbourlist = [[] for n in range(self.natoms)]
		for pair in pairs:
			i = pair[0]
			j = pair[1]
			self.neighbourlist[i].append(j)
			self.neighbourlist[j].append(i)
	def load_atoms(self):
		"""
		Reads `pdb_file`, sets `atoms` list and `natoms`. 
		Only heavy atoms are retained.
		"""
		with open(self.pdb_in) as f:
			pdb_lines = f.readlines()
		for line in pdb_lines:
			if line.startswith("ATOM"):
				atom = pdbatom.PdbAtom()
				atom.parse(line)
				if not (atom.elementSymbol=='H'):
					self.atoms.append(atom)
		self.natoms = len(self.atoms)
	def generate_sphere_points(self,n):
		"""
		Returns list of coordinates on a sphere of radius 1 using the Golden-Section Spiral algorithm.
		
		:param n: number of dots on the sphere
		:type n: `int`
		
		:return: list of coordinates of points covering a sphere of radius 1
		:rtype: `numpy.array` shape (`n`,3)
		
		"""
		points = numpy.zeros((n,3))
		inc = math.pi * (3 - math.sqrt(5))
		offset = 2 / float(n)
		for k in range(n):
			y = k * offset - 1 + (offset / 2)
			r = math.sqrt(1 - y*y)
			phi = k * inc
			points[k,:] = numpy.array([math.cos(phi)*r, y, math.sin(phi)*r])
		return points
	def compute_sas(self):
		"""
		Finds which dots on the atom spheres are solvent-acessible. Sets 
		`indices_selected_dots_per_atom`. Also sets the total number of 
		solvent accessible atom dots for the entire structure, this is
		used by `export_xyz`
		"""
		self.ndots_sas = 0
		for n in range(self.natoms):
			dots = self.atom_position[n,:]+ (self.atom_radius[n]+self.probe_radius)*self.generate_sphere_points(self.dots_per_atom)
			self.atom_dots.append(dots)
			exclude_dots = set()
			keep_dots = set([n for n in range(self.dots_per_atom)])
			kdtree = scipy.spatial.KDTree(dots)
			for nn in self.neighbourlist[n]:
				rnn = self.atom_position[nn]
				cutoff = self.atom_radius[nn]+self.probe_radius
				excl = set(kdtree.query_ball_point(rnn,cutoff))
				exclude_dots = exclude_dots.union(excl)
			keep_dots = keep_dots.difference(exclude_dots)
			self.indices_selected_dots_per_atom.append(list(keep_dots))
			self.ndots_sas+=len(keep_dots)
	def export_xyz(self,filename):
		"""
		Exports solvent accessible dots to xyz file `filename`
		
		:param filename: path to xyz file
		:type filename: `str`
		
		"""
		f = open(filename,'w')
		f.write(str(self.ndots_sas)+"\n")
		f.write("sas for "+self.pdb_in+"\n")
		for n in range(self.natoms):
			dot_indices = self.indices_selected_dots_per_atom[n]
			for i in dot_indices:
				dot = self.atom_dots[n][i,:]
				x = dot[0]
				y = dot[1]
				z = dot[2]
				f.write("dot   "+str(x)+"   "+str(y)+"   "+str(z)+"\n")
