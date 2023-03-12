""" 

module containing the class Dx

"""

import gridData 
import numpy
import scipy.interpolate

class Dx:
	"""Reads .dx grid data files produced by the Poisson-Boltzmann solver APBS 
	(using module gridData), and interpolates values for electrostatic potentials
	from the grid data. Dependencies: gridData, numpy, scipy.interpolate 
	
	:param dx_file: path to dx file to be read
	:type dx_file: `str`
	
	:param dxdata: datastructure containing the dx data
	:type dxdata: gridData.Grid
	
	:param value: interpolating function for the grid data, use as dx.value(numpy.array[x,y,z]).
	:type value: scipy.interpolate.RegularGridInterpolator
	
	:Example:
	
	.. code-block:: python3
	
		# lets get the electrostatic potential at some position r from a dx file produced by APBS
		dx = dx.Dx("pot.dx")
		r = numpy.array([3.2,4.5,4.6])
		print( "psi(r) = ", dx.value(r))	
	
	
"""
	def __init__(self,dx_file):
		"""
		Constructor. Reads the .dx file and sets up the interpolator `value`.
		
		:param dx_file: path to dx file
		:type dx_file: `str` 
		
		"""
		self.dx_file = dx_file
		self.dx = gridData.Grid(self.dx_file)
		self.dxdata = self.dx.grid
		self.value = scipy.interpolate.RegularGridInterpolator((self.dx.midpoints[0], self.dx.midpoints[1],self.dx.midpoints[2]), self.dxdata)




		
		

