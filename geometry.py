"""

convenience geometry functions for pyrosetta poses

"""

import math
import numpy
import pyrosetta

def get_rcm(pose):
	"""
	convenience function returning center of mass of pose as `numpy.array`
	
	:param pose: input pose
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
	
	:return: center of mass of pose
	:rtype: `numpy.array` shape (3,)
		
	"""
	nres = pyrosetta.rosetta.core.pose.nres_protein(pose)
	rcm_xyz = pyrosetta.rosetta.core.pose.center_of_mass(pose, 1,nres)
	return numpy.array([rcm_xyz.x,rcm_xyz.y,rcm_xyz.z])

def rotate_x(pose, theta_deg):
	"""
	rotates pose around x-axis over theta_deg degrees
	
	:param pose: input pose
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
		
	"""
	theta = (theta_deg/360.0)*2*math.pi
	Rx = numpy.matrix([[ 1, 0           , 0           ],
                   [ 0, math.cos(theta),-math.sin(theta)],
                   [ 0, math.sin(theta), math.cos(theta)]])
	pose.rotate(Rx)
  
  
def rotate_y(pose, theta_deg):
	"""
	rotates pose around y-axis over theta_deg degrees
	
	:param pose: input pose
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
		
	"""
	theta = (theta_deg/360.0)*2*math.pi
	Ry =  numpy.matrix([[ math.cos(theta), 0, math.sin(theta)],
				[ 0           , 1, 0           ],
                [-math.sin(theta), 0, math.cos(theta)]])      
	pose.rotate(Ry)

def rotate_z(pose,theta_deg):
	"""
	rotates pose around x-axis over theta_deg degrees
	
	:param pose: input pose
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
		
	"""
	theta = (theta_deg/360.0)*2*math.pi
	Rz = numpy.matrix([[ math.cos(theta), -math.sin(theta), 0 ],
                   [ math.sin(theta), math.cos(theta) , 0 ],
				   [ 0           , 0            , 1 ]])
	pose.rotate(Rz)
