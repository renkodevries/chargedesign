r"""
generates relevant residuesets for later use in script chargedesign.
Call as: 

.. code-block:: python3

	python makeresiduesets.py

Inputfile is `input/makeresiduesets.json`, sample imputfile:

.. code-block:: javascript
	
	{
	"pdb_symm" : "pdb/in/M_C3.pdb",
	"contact_distance_interface" : 7,
	"residuesets_out_filename" : "residuesets/M_residuesets.json",
	"pml_out_filename" : "pml/M_designable.pml",
	"composition_out_filename" : "residuesets/M_comp.json",
	"exclude_residues" : [1,2,3,4,5,6,273,274,275,276,277,278,279]
	}

Where:

* `pdb_symm` is  a pdb file containing the structure of a Cn symmetric multimer
* `contact_distance_interface` is the contact distance in Angstroms used to determine the interface between chains A and B of the multimer
* `residuesets_out_filename` is the filename for the output residueset .json file
* `pml_out_filename` is the filename for the output .pml PyMol scriptfile for selecting the residuesets (run from within PyMol), 
* `composition_out_filename` is the filename for the .json outputfile containing info on the composition of the output residuesets. 
* `exclude_residues` is a residueset to be excluded from the final `designable` residueset produced by the script. 

Note that residuenumbers are the pyrosetta residuenumbers for the monomer.
	
	
"""
import pyrosetta
import sys
import json
import os

# local modules
import residuesets

# main
if __name__ == "__main__":
	scriptname = os.path.basename(__file__).rsplit(".",maxsplit=1)[0]
	thisscript = "gets residuesets for subsequent surface design"
	inputdir = "input"
	inputfile = os.path.join(inputdir,scriptname+".json")
	
	print("-------------------------")
	print("script:", scriptname)
	print(thisscript)
	print("inputfile: ",inputfile)
	print("-------------------------")

	print('\n ====> read & parse json input file\n')
	with open(inputfile, 'r') as f: input = json.load(f)
	pdb_symm = input["pdb_symm"]
	contact_distance_interface = input["contact_distance_interface"]
	residuesets_out_filename = input["residuesets_out_filename"]
	pml_out_filename = input["pml_out_filename"]
	composition_out_filename = input["composition_out_filename"]
	exclude_residues = input["exclude_residues"]
	
	print("\n ====> initialize pyRosetta, load symmetric input pdb\n ")
	print("...pdb_symm = ",pdb_symm)

	pyrosetta.init()
	pose_symm = pyrosetta.io.pose_from_pdb(pdb_symm)

	print("\n ====> make monomer (chain A) and dimer (chain A,B) poses, save as pdb\n ")
	pose_monomer =  pyrosetta.rosetta.core.pose.Pose()
	nresidues_monomer = pose_monomer.total_residue()
	pose_dimer = pyrosetta.rosetta.core.pose.Pose()
	chain_A_id = pyrosetta.rosetta.core.pose.get_chain_id_from_chain('A',pose_symm)
	chain_B_id = pyrosetta.rosetta.core.pose.get_chain_id_from_chain('B',pose_symm)
	pyrosetta.rosetta.core.pose.append_pose_to_pose(pose_monomer, pose_symm.split_by_chain(chain_A_id))
	pyrosetta.rosetta.core.pose.append_pose_to_pose(pose_dimer, pose_symm.split_by_chain(chain_A_id))
	pyrosetta.rosetta.core.pose.append_pose_to_pose(pose_dimer, pose_symm.split_by_chain(chain_B_id))
	pdb_no_ext = pdb_symm.rsplit(".",maxsplit=1)[0]
	pdb_monomer = pdb_no_ext+"_A.pdb"
	pose_monomer.dump_pdb(pdb_monomer)
	pose_monomer = pyrosetta.io.pose_from_pdb(pdb_monomer)
	pdb_dimer = pdb_no_ext+"_A_B.pdb"
	pose_dimer.dump_pdb(pdb_dimer)
	pose_dimer = pyrosetta.io.pose_from_pdb(pdb_dimer)

	print("\n ====> get residueset surface\n")
	ressets = dict()
	ressets['all_surface']  = residuesets.residueset_surface(pose_monomer)

	print("\n ====> get residuesets interface_1 and interface_2\n")
	threshold_distance =  input["contact_distance_interface"]
	interface = residuesets.residueset_interface(pose_dimer, 1,2, threshold_distance)
	ressets['interface_1']  = interface["contacts_on_1"]
	ressets['interface_2']  = interface["contacts_on_2"]
	nresidues_monomer = pose_monomer.total_residue()
	ressets['interface_2']  = [n - nresidues_monomer for n in ressets['interface_2']]

	print("\n ====> get final residueset for surface design\n")
	ressets["all"] = [n+1 for n in range(nresidues_monomer)]
	ressets["exclude_residues"] = exclude_residues
	designable = set(ressets['all_surface'])
	designable = designable.difference(set(ressets['interface_1']))
	designable = designable.difference(set(ressets['interface_2']))
	designable = designable.difference(set(ressets['exclude_residues']))
	ressets["designable"] = list(designable)

	print("\n ====> save residuesets, pymol selections and residuesets compositions\n")
	print("...residuesets: ")
	for key in ressets.keys():
		print("	", key)

	# residuesets
	with open(residuesets_out_filename, 'w') as outfile:
		json.dump(ressets, outfile, indent=2)

	# pymol selections
	with open(pml_out_filename,'w') as outfile:
		for key in ressets:
			pymol_line = residuesets.residueset_to_pymolselection(key,ressets[key])
			outfile.write(pymol_line.strip()+"\n")

	# composition
	composition = dict()
	for key in ressets:
		composition[key] = residuesets.residueset_to_composition(pose_monomer,ressets[key])
	with open(composition_out_filename,'w') as outfile:
		json.dump(composition, outfile, indent=2)
	
	print("\n...residuesets_out_filename = ",residuesets_out_filename)
	print("...pml_out_filename = ",pml_out_filename)
	print("...composition_out_filename = ",composition_out_filename)
	print("done.")

