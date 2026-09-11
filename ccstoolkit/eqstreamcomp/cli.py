#!/usr/bin/python3

#===================================================================================================================
#---------------------------------------------------------------------------------------Import
#===================================================================================================================
import argparse
from tabulate import tabulate
import matplotlib.pyplot as plt

from .composition import get_composition
from .stoichiometry_map import get_stoichiometry_map
from .stability_map import get_stability_map
from ._reactions import get_reactions
from ccstoolkit.common._substances import get_substances_TD_data

#===================================================================================================================
#---------------------------------------------------------------------------------------Methods
#===================================================================================================================
def _print_output(output, output_file=None):
	if output_file:
		with open(output_file, "w") as f:
			f.write(output)
		print(f"Output written to {output_file}")
	else:
		print(output)

#===================================================================================================================
#---------------------------------------------------------------------------------------Main
#===================================================================================================================
def main():
	#---------------------------------------------------------------------------------------Parse the input
	parser = argparse.ArgumentParser(description="Compute the equilibrium composition.")
	
	#Print TD constants
	parser.add_argument("-c","--constants",action="store_true",help="Print the TD constants.")
	
	#Print reaction list
	parser.add_argument("-r","--reactions",action="store_true",help="Print the list of reactions.")
	
	#Calculate the equilibrium composition
	parser.add_argument(
		"-eq",
		action="store",
		type=float,
		nargs=6,
		metavar=("H", "N", "O", "S", "CO2", "T"),
		help="Calculate the equilibrium composition at a given H content [mol/m^3], N content [mol/m^3], O content [mol/m^3], S content [mol/m^3], CO2 activity [mol/m^3], and temperature [K]."
	)
	
	#Calculate the stoichiometry map
	parser.add_argument("-sto",action="store",type=float,metavar="ratio",help="Calculate the stoichiometry map at a given N/S ratio.")
	
	#Calculate the stability map
	parser.add_argument(
		"-sta",
		action="store",
		type=float,
		nargs=4,
		metavar=("N", "S", "CO2", "T"),
		help="Calculate the stability map at a given N content [mol/m^3], S content [mol/m^3], CO2 activity [mol/m^3], and temperature [K]."
	)
	
	#Print data
	parser.add_argument("-pr","--print",action="store_true",help="Print data.")
	
	#Output
	parser.add_argument("-o",action="store",type=str,metavar='file',help="Save text output to file.")
	
	
	args = parser.parse_args()
	
	#---------------------------------------------------------------------------------------Calculate
	#Calculate the equilibrium composition	
	if args.eq:
		H, N, O, S, CO2, T = args.eq
		c0 = {'H': H, 'N': N, 'O': O, 'S': S, 'CO2': CO2, 'T': T}
		
		c = get_composition(c0)
	
	#Calculate the stoichiometry map
	if args.sto:
		stoichiometry = get_stoichiometry_map({'N/S': args.sto})
		
	#Calculate the stability map	
	if args.sta:
		N, S, CO2, T = args.sta
		P = {'N': N, 'S': S, 'CO2': CO2, 'T': T}
		
		stability = get_stability_map(P)
		
	#---------------------------------------------------------------------------------------Output
	output = ''
	
	#Print TD constants
	if args.constants:
		substances = get_substances_TD_data()
		
		output += f"----------------------------------TD data-----------------------------------" + '\n'
		
		headers = ["Substance", "dfg [kJ/mol]", "dfh [kJ/mol]", "cp [J/mol/K]", "solid"]
		rows = [[key, data['dfg'], data['dfh'], data['cp'], data['solid']] for key,data in substances.items()]
		table = tabulate(rows, headers=headers, tablefmt="github", floatfmt=".2f") #tablefmt="grid"
		
		output += table + '\n\n'
	
	#Print reaction list
	if args.reactions:
		reactions = get_reactions()
		
		output += f"---------------------------------Reactions----------------------------------" + '\n'
		for key, reaction in reactions.items():
			output += reaction['reaction']['reaction'] + '\n'
			output += 'key: ' + key + '\n'
			output += 'drg: ' + f"{reaction['drg']*1e-3:.2f} kJ/mol" + '\n'
			output += 'drh: ' + f"{reaction['drh']*1e-3:.2f} kJ/mol" + '\n'
			output += 'drcp: ' + f"{reaction['drcp']:.2f} J/mol/K" + '\n'
			output += 'K_chi_298.15: ' + f"{reaction['K_chi_298']:.2e}" + '\n\n'
			
	#Equilibrium data
	if (args.print or args.o):
	
		if args.eq and type(c)!=int:
			output += f"--------------------------Equilibrium composition---------------------------" + '\n'
			
			headers = ["Substance", "c [mol/m^3]"]
			rows = [[key, value] for key, value in c.items()]
			table = tabulate(rows, headers=headers, tablefmt="github", floatfmt=".2e")
			
			#for key, value in c.items():
			#	output += key + ":\t" + f"{value:.2e} mol/m^3\n"
				
			output += table+'\n\n'
				
		if args.sto:
			output += f"-----------------------------Stoichiometry map------------------------------" + '\n'
			for region in stoichiometry:
				output += region['name'] + '\n'
				output += f"Area: {region['area']:.2f}" + '\n'
				output += f"Centroid: {list(region['centroid'])}" + '\n'
				output += f"Vertices: {region['points']}" + '\n\n'
				
		if args.sta:
			output += f"-------------------------------Stability map--------------------------------" + '\n'
			for region in stability:
				output += region['name'] + '\n'
				output += f"Area: {region['area']:.2f}" + '\n'
				output += f"Centroid: {list(region['centroid'])}" + '\n'
				output += f"Vertices: {region['points']}" + '\n\n'
			
	#Print
	if output:
		_print_output(output, args.o)
	
	return 0

if __name__ == "__main__":
	main()
