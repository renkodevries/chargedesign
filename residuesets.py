""" 

module with functions for selecting residuesets and functions operating on residuesets,
where residuesets are lists of pyrosetta residuenumbers. 

"""

import pyrosetta
import math
import numpy

def residueset_to_composition(pose,residueset):
	"""
	given a pose and residueset, returns dictionary with 
	information about the aminoacid composition of the residueset
	
	
	:param pose: input pose 
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
	
	:param residueset: input residueset (list of pyrosetta residuenumbers) 
	:type residueset: `list` of `int`
	
	:return: dictionary with information about the aminoacid composition of the residueset
	:rtype: `dict`
	
	"""
	sequence = pose.sequence()
	composition = dict()
	charge = 0
	for n in range(len(sequence)):
		if n in residueset:
			aa = sequence[n]
			# charge
			if (aa=='K' or aa=='R'): charge+=1
			elif (aa=='E') or (aa=='D'): charge-=1
			# composition
			if aa in composition:
				composition[aa]+=1
			else:
				composition[aa]=1
	data = dict()
	data['residueset_length'] = len(residueset)
	data['composition'] = composition
	data['charge'] = charge
	return data

def residueset_to_charge(residueset,pose):
	"""
	given a pose and residueset, returns net charge of residueset
	
	:param pose: input pose 
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
	
	:param residueset: input residueset (list of pyrosetta residuenumbers) 
	:type residueset: `list` of `int`
	
	:return: net charge of residueset
	:rtype: `int`
	
	"""
	sequence = pose.sequence()
	Q = 0
	for i in residueset:
		aa = sequence[i]
		if aa in ['D','E']: Q-=1
		elif aa in ['K','R']: Q+=1
	return Q

def residueset_to_resfile(filename,residueset_nataa, residueset_allaa):
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
					
def residueset_to_pyrosettaselector(residueset):
	"""
	converts residueset to rosetta residue_selector
	
	:param residueset: residueset to be converted 
	:type residueset: `list` of `int`
	
	:return: residue_selector
	:rtype: `pyrosetta.rosetta.utility.vector1_bool`
	"""
	sep=","
	selection_str = sep.join([str(i) for i in residueset])
	selector =  pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector(selection_str) 
	return selector

def pyrosettaselector_to_residueset(selector,ntot):
	
	"""
	
	converts rosetta residue_selector to residueset  
	
	:param selector: residue_selector to be converted 
	:type selector: `pyrosetta.rosetta.utility.vector1_bool`
	
	:param ntot: total number of residues of the residue_selector - NOT CLEAR WHY THIS IS NEEDED
	:type ntot: `int`

	:return: residueset
	:rtype: `list` of `int`
	
	"""
	
	res_set = set()
	for n in range(ntot):
		if selector[n+1]: res_set.add(n+1)
	return list(res_set)

def residueset_to_pyrosettamask(residueset,nres):
	"""
	converts residueset to rosetta residue_selector (or mask) - SEEMS TO BE DUPLICATE OF residueset_to_pyrosettaselector 

	:param residueset: residueset to be converted 
	:type residueset: `list` of `int`

	:param nres: total number of residues of the residue_selector 
	:type nres: `int`

	:return: residue_selector
	:rtype: `pyrosetta.rosetta.utility.vector1_bool`
	"""
	mask = pyrosetta.rosetta.utility.vector1_bool()
	for i in range(nres):
		ires = i+1
		if ires in residueset:
			mask.append(True)
		else:
			mask.append(False)
	return mask

def residueset_close_to_target_residueset(pose,target_residueset,distance_threshold):
	
	"""
	
	find residueset of residues close to residues in a target_residueset given some distance_treshold 

	:param pose: input pose 
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
	
	:param target_residueset: residueset of residues close to which we are seeking other residues 
	:type target_residueset: `list` of `int`

	:param distance_treshold: distance_treshold in Angstrom 
	:type distance_treshold: `float`

	:return: residueset
	:rtype: `list` of `int`
	
	"""	
	
	nresidues = pose.total_residue() 
	index_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
	index_string = ",".join([str(n) for n in target_residueset])
	index_selector.set_index(index_string)
	contact_selector = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
	contact_selector.threshold(distance_threshold)
	contact_selector.central_residue_group_selector(index_selector)
	contacts_residueset = pyrosettaselector_to_residueset(contact_selector.apply(pose),nresidues)
	return contacts_residueset

def residueset_to_pymolselection(selectionname, residueset):
	"""
	convert residueset to pymol scriptline for selecting these residues

	:param selectionname: name for pymol selection 
	:type selectionname: `str`
	
	:param residueset: residueset to be converted 
	:type residueset: `list` of `int`

	:return: pymol script line selecting the residueset and naming it `selectionname`
	:rtype: `str`

	"""	
	line="select "+selectionname+","
	for n in range(len(residueset)):
		resno = residueset[n]
		if (n==0):
			line+=" resi "+str(resno)
		else:
			line+=" + resi "+str(resno)
	return line

def residueset_surface(pose):
	"""
	get residueset of surface residues for given pose
	
	:param pose: input pose 
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
	
	:return: residueset of surface residues
	:rtype: `list` of `int`
	
	"""	
	nresidues = pose.total_residue() 
	layer_selector = pyrosetta.rosetta.core.select.residue_selector.LayerSelector()
	layer_selector.set_layers(False, False, True)
	surface_selector = layer_selector.apply(pose)
	return  pyrosettaselector_to_residueset(surface_selector,nresidues)

def residueset_interface(pose, chain1_id, chain2_id, threshold_distance):
	"""
	get residueset of residues forming the interface between chain_1 and chain_2 of given pose
	
	:param pose: input pose 
	:type pose: `pyrosetta.rosetta.core.pose.Pose`
	
	:param chain1_id: pyrosetta chain number of chain_1 
	:type pose: `int`
	
	:param chain2_id: pyrosetta chain number of chain_2 
	:type pose: `int`
	
	:return: dictionary with as keys `contacts_on_1` and `contacts_on_2` 
	and as values the residuesets forming the interface between chain_1 and chain_2
	:rtype: `dict`
	
	"""	
	nresidues = pose.total_residue()
	chain_1 = pyrosetta.rosetta.core.pose.get_chain_from_chain_id(1,pose)
	chain_2 = pyrosetta.rosetta.core.pose.get_chain_from_chain_id(2,pose)
	chain_1_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain_1)
	chain_2_selector = pyrosetta.rosetta.core.select.residue_selector.ChainSelector(chain_2)
	chain_1_set = set(pyrosettaselector_to_residueset(chain_1_selector.apply(pose), nresidues))
	chain_2_set = set(pyrosettaselector_to_residueset(chain_2_selector.apply(pose), nresidues))
	contact_selector = pyrosetta.rosetta.core.select.residue_selector.CloseContactResidueSelector()
	contact_selector.threshold(threshold_distance)
	contact_selector.central_residue_group_selector(chain_1_selector)
	contacts_1_2_set = set(pyrosettaselector_to_residueset(contact_selector.apply(pose),nresidues))
	contact_selector.central_residue_group_selector(chain_2_selector)
	contacts_2_1_set = set(pyrosettaselector_to_residueset(contact_selector.apply(pose), nresidues))
	contacts_on_2 = list(contacts_1_2_set.difference(chain_1_set))
	contacts_on_1 = list(contacts_2_1_set.difference(chain_2_set))
	out = dict()
	out["contacts_on_1"]= list(contacts_on_1)
	out["contacts_on_2"]= list(contacts_on_2)
	return out
