#!/usr/bin/python3

#===================================================================================================================
#---------------------------------------------------------------------------------------Import
#===================================================================================================================
import re
import sys
import argparse
from tabulate import tabulate
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

from .maps import get_maps
from ._reactions import get_domain
from ._reactions import get_bounds
from ._reactions import get_reactions
from ._substances import get_substances_TD_data

#===================================================================================================================
#---------------------------------------------------------------------------------------Constants
#===================================================================================================================
_corr_products_colors = {
    'Fe': 'tab:blue',
    'FeCO3': 'tab:orange',
    'Fe(OH)2': 'navy',
    'FeS2': 'tab:green',
    'FeS': 'tab:orange',
    'FeSO4': 'tab:red',
    'FeSO4.H2O': 'tab:purple',
    'FeSO4.7H2O': 'tab:pink',
    'Fe(NO3)2': 'tab:olive',
    'Fe3O4': 'tab:green',
    'Fe2O3': 'tab:red',
    'FeO(OH)': 'tab:purple',
    'Fe2(SO4)3': 'tab:olive'
}

#===================================================================================================================
#---------------------------------------------------------------------------------------Methods
#===================================================================================================================
def print_output(output, output_file=None):
	if output_file:
		with open(output_file, "w") as f:
			f.write(output)
		print(f"Output written to {output_file}")
	else:
		print(output)
		
def verify_domain(P):
	domain = get_domain()
		
	if not domain['S']['min'] <= P['S'] <= domain['S']['max']:
		print(f"S out of bounds; [{domain['S']['min']:.3f} mM : {domain['S']['max']:.3f} mM]")
		sys.exit(1)
			
	if not domain['N']['min'] <= P['N'] <= domain['N']['max']:
		print(f"N out of bounds; [{domain['N']['min']:.3f} mM : {domain['N']['max']:.3f} mM]")
		sys.exit(1)
			
	if not domain['C']['min'] <= P['C'] <= domain['C']['max']:
		print(f"CO2 out of bounds; [{domain['C']['min']:.3f} mM : {domain['C']['max']:.3f} mM]")
		sys.exit(1)
			
	if not domain['T']['min'] <= P['T'] <= domain['T']['max']:
		print(f"T out of bounds; [{domain['T']['min']:.3f} K : {domain['T']['max']:.3f} K]")
		sys.exit(1)
		
	return True
	
def plot(maps, fname=None):
	Bounds = get_bounds()

	fig, axs = plt.subplot_mosaic(
		[["O", "C"],["N", "S"],["legend", "legend"]],
		figsize=(9, 8),
		height_ratios=[1, 1, 0.3],
		layout="tight"
	)
	
	for ax, (key,map_) in zip(axs.values(),maps.items()):
		for region in map_:
			label = re.sub(r'(?<!\.)\d', r'$_\g<0>$', region['name'])
			ax.fill(*zip(*region['points']),hatch='xxxx',color='white',edgecolor=_corr_products_colors[region['name']],label=label)
		
		ax.set_xlim(*Bounds['x'])
		ax.set_ylim(*Bounds['y'])
	
		ax.set_title(f"{key} map",fontsize=16)
		
		ax.set_xlabel(r'$\lg\text{H}_2\text{O}$', fontsize=12)
		ax.set_ylabel(r'$\lg\text{O}_2$', fontsize=12)
		
		
	axs['legend'].axis("off")
	
	handles = {k: j for i in ('O','C','N','S') for j,k in zip(*axs[i].get_legend_handles_labels())}
	labels = set([j for i in ('O','C','N','S') for j in axs[i].get_legend_handles_labels()[1]])
	axs["legend"].legend([handles[label] for label in labels], labels, loc="center", ncol=5, frameon=False, fontsize=12)
	
	if fname:
		plt.savefig(fname+'.png', format='png') #, transparent=True
		print(f"Plot saved to {fname}.png")
	else:
		plt.show()

#===================================================================================================================
#---------------------------------------------------------------------------------------Main
#===================================================================================================================
def main():
	#---------------------------------------------------------------------------------------Parse the input
	parser = argparse.ArgumentParser(description="Compute corrosion products maps.")
	
	#Print TD constants
	parser.add_argument("-c","--constants",action="store_true",help="Print the TD constants.")
	
	#Print reaction list
	parser.add_argument("-r","--reactions",action="store_true",help="Print the list of reactions.")
	
	#Calculate the phase diagrams
	parser.add_argument(
		"-p",
		action="store",
		type=float,
		nargs=4,
		metavar=("S", "N", "C", "T"),
		help="Calculate the phase diagrams at a given total S [mM], total N [mM], C activity [mM], temperature [K]."
	)
	
	#Print reaction list
	parser.add_argument("-pr","--print",action="store_true",help="Print the maps.")
	
	#Plot
	parser.add_argument("-pl","--plot",action="store_true",help="Plot the maps.")
	
	#Save plot
	parser.add_argument("-sp",action="store",type=str,metavar='file',help="Save the plot.")
	
	#Output
	parser.add_argument("-o",action="store",type=str,metavar='file',help="Output file name.")
	
	args = parser.parse_args()
	
	#---------------------------------------------------------------------------------------Calculate
	#Calculate the phase diagrams	
	if args.p:
		S, N, C, T = args.p
		P = {'S': S, 'N': N, 'C': C, 'T': T}
		
		if verify_domain(P):
			maps = get_maps(P)
	
	#---------------------------------------------------------------------------------------Output
	output = ''
	
	#Print TD constants
	if args.constants:
		substances = get_substances_TD_data()
		
		headers = ["Substance", "dfg [kJ/mol]", "dfh [kJ/mol]", "cp [J/mol/K]", "solid"]
		rows = [[key, data['dfg'], data['dfh'], data['cp'], data['solid']] for key,data in substances.items()]
		table = tabulate(rows, headers=headers, tablefmt="github", floatfmt=".2f") #tablefmt="grid"
		
		output += table + '\n'
	
	#Print reaction list
	if args.reactions:
		reactions = get_reactions()
		
		output += '\n'
		for key, reaction in reactions.items():
			output += reaction['reaction']['reaction'] + '\n'
			output += 'key: ' + key + '\n'
			output += 'drg: ' + f"{reaction['drg']*1e-3:.2f} kJ/mol" + '\n'
			output += 'drh: ' + f"{reaction['drh']*1e-3:.2f} kJ/mol" + '\n'
			output += 'drcp: ' + f"{reaction['drcp']:.2f} J/mol/K" + '\n'
			output += 'K_chi_298.15: ' + f"{reaction['K_chi_298']:.2e}" + '\n\n'
	
	#Print reaction list
	if (args.print or args.o) and args.p:
		output += '\n'
		for key, map_ in maps.items():
			output += f"------------------------------------{key} map------------------------------------" + '\n'
			for region in map_:
				output += region['name'] + '\n'
				output += f"Area: {region['area']:.2f}" + '\n'
				output += f"Centroid: {region['centroid']}" + '\n'
				output += f"Vertices: {region['points']}" + '\n\n'
	
	#Print
	if output:
		print_output(output, args.o)
	
	#Print reaction list
	if (args.plot or args.sp) and args.p:
		plot(maps,args.sp)
	
	return 0

if __name__ == "__main__":
	main()
