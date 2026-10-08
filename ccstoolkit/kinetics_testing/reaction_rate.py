
## functions for calculating rates of elementary and net reactions

import numpy as np
import re


# rate of an elementary reaction
def elementary_reaction_rate(c,coef,k):
	# c = list of concentrations e.g. [1,2,0.5,...]
	# coef = list of stoichiometric coefficients e.g. [0,0,-1,1,...]
	# k = rate constant e.g. 1,2e-5
		# rate = k * c[1]**(-coef[1]) * c[2]**(-coef[2]) * ... for coef[i]<0

	# create an array of 1's
	rate = np.ones_like(c, dtype=float)

	# determine where the stoichiometric coefficient is negative, i.e reactant
	mask = coef < 0

	# raise the concentration to the power of the stoichiometric coefficient
	rate[mask] = c[mask]**(-coef[mask])

	# multiply by the rate constant
	rate = k*np.prod(rate)

	return rate

# rate of a net reaction
def net_reaction_rate(c, rate_law, k, species_list):
	# c = list of concentrations e.g. [1,2,0.5,...]
	# parsed rate = 
	# k = rate constant e.g. 1,2e-5

	if rate_law == "k":
		rate = k
		return rate

	# remove whitespace
	rate_law = rate_law.strip()
	rate_law = rate_law.replace(" ", "")

	# check rate starts with k*, and then removes it
	if not rate_law.startswith("k*"):
		raise ValueError(f"Bad rate format '{rate_law}'. Rate must start with 'k*', e.g. k*[A]*[B]^2. If the rate is only dependent on the rate constant please input it as k. If the reaction is reversible please input it as 'k' in both rate laws.")
	rate_concs = rate_law[2:] # gives e.g. [A]*[B]^2

	# creates an array of 1s
	rate = np.ones_like(c, dtype=float)

	# splits expression via '['
	concs = rate_concs.split("[")

	# gets species and coefficients
	for conc in concs[1:]:
		# check if corresponding ']' is present
		if "]" not in conc:
			raise ValueError(f"Missing ']' in rate expression '{rate_law}'.")

		species, coef_str = conc.split("]") # species = e.g. 'A', coef_str = e.g. '^2'

		# check species is in species_list
		if species not in [spec['name'] for spec in species_list]:
			raise ValueError(f"Species '{species}' in rate expression '{rate_law}' is not given as a species in the input data. Given species: {', '.join(spec['name'] for spec in species_list)} ")

		# Find the coefficient
		coef_str = coef_str.rstrip("*") 
		match = re.fullmatch(r"\^(\d*\.?\d+(?:/\d*\.?\d+)?)", coef_str)
		if coef_str == "":
			power = 1.0
		elif match:
			coef_str = match.group(1)
			if "/" in coef_str:
				num, den = coef_str.split("/")
				power = float(num) / float(den)
			else:
				power = float(coef_str)
		else:
			raise ValueError(f"Could not parse power in rate expression '{rate_law}'. Please use integers, decimals or fractions, e.g. ^2, ^0.5, ^3/2.")

		# find position of species in the list
		index = None
		for i, spec in enumerate(species_list):
			if spec["name"] == species:
				index = i
				break

		# calculate concentrations
		rate[index] = c[index]**power

	rate = k * np.prod(rate)

	return rate

