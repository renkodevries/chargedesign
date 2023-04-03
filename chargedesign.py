import pyrosetta
import json
import math
import string
import os
import sys
import numpy
import random
import argparse

# local
import potentials
import chargedesigner

if __name__ == "__main__":

	print ("\n ====> initialize pyRosetta \n")
	pyrosetta.init()

	print('\n ====> read json input file\n')
	thisscript = "chargedesign: find surface residues that minimize fluctuations of surface electrostatic potential."
	parser = argparse.ArgumentParser(description = thisscript)
	parser.add_argument('json_in',  help = "json input file")
	args = parser.parse_args()
		
	print("\n ===>  initialize instance of ChargeDesigner\n")
	print("... input file: ",args.json_in)
	chargedesigner = chargedesigner.ChargeDesigner(args.json_in)
	
	print("\n ===>  do simulated annealing steps\n")
	print("...",chargedesigner.n_steps_hold_start," steps at kT_sampling = ",chargedesigner.kT_sampling_start)
	print("...",chargedesigner.n_steps_linear_gradient," steps, linear gradient from kT_sampling = ",chargedesigner.kT_sampling_start, "...",chargedesigner.kT_sampling_end)
	print("...",chargedesigner.n_steps_hold_end," steps at kT_sampling = ",chargedesigner.kT_sampling_end)
	
	while chargedesigner.step < chargedesigner.n_total_steps:
		
		chargedesigner.step+=1
		print("\n ===============")
		print("    step = ",chargedesigner.step)
		print(" ===============\n")	
		
		print("\n ##### set kT_sampling\n")
		chargedesigner.set_kT_sampling()
		
		print("\n ##### do trial move\n")
		chargedesigner.do_trial_move()
		
		print("\n ##### calculate trial move scores\n")
		chargedesigner.set_scores()
		
		print("\n ##### apply acceptance criterium\n")
		chargedesigner.apply_acceptance_criterium()
		
		print("\n ##### save scores\n")
		chargedesigner.save_scores()
		
		print("\n ##### check for new minimum\n")
		chargedesigner.check_for_new_minimum()
		
	print("done with chargedesign.") 


