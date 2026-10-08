
## PARSE_REACTION parses a reaction based on 3 conditions: 
	## irreversible/ reversible
	## elementary/ net 
	## rate constants based on k or A, Ea

import numpy as np
import re

class PARSE_REACTION():

	def __init__(self, line, reaction_type, parameter_type, species, temperature, calculate_k_function):
		self._line = line
		self._species = species
		self._species_index = {s["name"]: i for i, s in enumerate(self._species)}
		self._temperature = float(temperature)
		self._calculate_k = calculate_k_function

		# checks reversibility of reaction based on arrow format
		reversibility = self._check_reversibility()

		# parses each reaction based on the 3 conditions
		if reversibility == "irreversible":
			if reaction_type == "elementary":
				if parameter_type == "rate constant":
					self._reaction = self._parse_irrev_elem_k()

				elif parameter_type == "arrhenius":
					self._reaction = self._parse_irrev_elem_AEa()

			elif reaction_type == "net":
				if parameter_type == "rate constant":
					self._reaction = self._parse_irrev_net_k()

				elif parameter_type == "arrhenius":
					self._reaction = self._parse_irrev_net_AEa()

		elif reversibility == "reversible":
			if reaction_type == "elementary":
				if parameter_type == "rate constant":
					self._reaction = self._parse_rev_elem_k()

				elif parameter_type == "arrhenius":
					self._reaction = self._parse_rev_elem_AEa()

			elif reaction_type == "net":
				if parameter_type == "rate constant":
					self._reaction = self._parse_rev_net_k()

				elif parameter_type == "arrhenius":
					self._reaction = self._parse_rev_net_AEa()

	# checks reversibility of reaction based on arrow format
	def _check_reversibility(self):
		if "->" in self._line:
			reversibility = "irreversible"
		elif "<=>" in self._line:
			reversibility = "reversible"
		elif "<->" in self._line:
			raise ValueError(f"'<->' is not a permitted reaction arrow. Please use '->' for an irreversible reaction or '<=>' for a reversible reaction.")
		else:
			error = "\n\nEnsure that the reaction arrow used is '->' or '<=>' depending on the reversibility of the reaction."
			raise ValueError(self._bad_reaction_format(error))

		return reversibility

	# error message for a bad reaction line format
	def _bad_reaction_format(self, error=""):
		return (
			f"Bad reaction line format in line: '{self._line}'."
			"\n\nElementary reactions are formatted as: 'reaction name | kinetic parameter(s) | reaction (e.g. A -> B or A <=> B)'."
			"\n\nNet reactions are formatted as: 'reaction name | rate law | kinetic parameters(s) | reaction (e.g. A -> B or A <=> B)'."
			f"{error}"
		)

	# gets stoichiometric coefiicients form equation
	def _find_coefficients(self, eq):
		lhs, rhs = eq.split("->")
		coeff = np.zeros(len(self._species))
		self._parse_side(lhs, -1, coeff)
		self._parse_side(rhs, +1, coeff)

		return coeff

	def _parse_side(self, side, sign, coeff):
		if not side.strip():
			raise ValueError("Empty reaction side detected.")
			
		for term in side.split("+"):
			match = re.match(
				r'^(\d+(?:/\d+|\.\d+)?)?\s*([A-Za-z].*)',
				term.strip()
			)

			if not match:
				raise ValueError(f"Cannot parse reaction term: '{term}'")

			coef_str, species = match.groups()

			if coef_str is None:
				coef_val = 1.0
			elif "/" in coef_str:
				num, den = coef_str.split("/")
				coef_val = float(num) / float(den)
			else:
				coef_val = float(coef_str)

			species = species.strip()

			if species not in self._species_index:
				raise ValueError(f"Unknown species '{species}' in equation. Please add this species to the input.")

			i = self._species_index[species]
			coeff[i] += sign * coef_val

	# counts the number of reactants in a reaction equation
	def _count_reactants(self, eq):
		lhs = eq.split("->")[0]
		total_reactants = 0.0

		for term in lhs.split("+"):
			match = re.match(r'^(\d+(?:/\d+|\.\d+)?)?\s*([A-Za-z].*)', term.strip())
			if not match:
				raise ValueError(f"Cannot parse reaction term: '{term}'")

			coef_str = match.group(1)
			if coef_str is None:
				total_reactants += 1.0
			elif "/" in coef_str:
				num, den = coef.str.split("/")
				total_reactants += float(num)/float(den)
			else:
				total_reactants += float(coef_str)

		return total_reactants

	# reaction = irreversible, elementary, kinetics = rate constant k
	def _parse_irrev_elem_k(self):
		# format of line: reaction name | rate constant | reaction
		reaction = []

		parts = self._line.split("|")
		if len(parts) != 3:
			error = (
				"\n\nFor an irreversible elementary reaction controlled by rate constant k, kinetic parameter(s) must be the rate constant (k) of the reaction."
				"\ne.g. r1 | 1e3 | A -> B"
			)
			raise ValueError(self._bad_reaction_format(error))

		name, k, eq = parts

		if float(k) <= 0:
			raise ValueError(f"Rate constant for reaction in line '{self._line}' must be positive")

		reac_parts = {
			"name": name.strip(),
			"rate constant": float(k),
			"equation": eq.strip(),
			"coefficients": self._find_coefficients(eq.strip()),
			"num reactants": self._count_reactants(eq.strip()),
		}

		reaction.append(reac_parts)
		return reaction 

	# reaction = irreversible, elementary, kinetics = arrhenius parameters A and Ea
	def _parse_irrev_elem_AEa(self):
		# format of line: reaction name | arrhenius parameter A | arrhenius parameter Ea | reaction
		reaction = []

		parts = self._line.split("|")
		if len(parts) != 4:
			error = (
				"\n\nFor an irreversible elementary reaction controlled by arrhenius parameters A and Ea, the kinetic parameter(s) must be A and Ea (in J /mol)."
				"\ne.g. r1 | 10e6 | 40000 | A -> B"
			)
			raise ValueError(self._bad_reaction_format(error))

		name, A, Ea, eq = parts

		if float(A) <= 0 or float(Ea) <= 0:
			raise ValueError(f"Arrhenius parameters A and Ea in line '{self._line}' must both be positive")

		reac_parts = {
			"name": name.strip(),
			"rate constant": self._calculate_k(float(A), float(Ea), self._temperature),
			"equation": eq.strip(),
			"coefficients": self._find_coefficients(eq.strip()),
			"A": A,
			"Ea": Ea,
			"num reactants": self._count_reactants(eq.strip()),
		}

		reaction.append(reac_parts)
		return reaction

	# reaction = irreversible, net, kinetics = rate constant k
	def _parse_irrev_net_k(self):
		# format of line: reaction name | rate law | rate constant | reaction
		reaction = []

		parts = self._line.split("|")
		if len(parts) != 4:
			error = (
				"\n\nFor an irreversible net reaction controlled by rate constant k, make sure the reaction line is formatted as above."
				"\nThe rate law must be of the format: k*[A]^2, where k is the rate constant and A (and additional species) are species concentrations."
				"\ne.g. r1 | k*[A]^2 | 1e3 | A -> B"
			)
			raise ValueError(self._bad_reaction_format(error))

		name, rate, k, eq = parts

		if float(k) <= 0:
			raise ValueError(f"Rate constant for reaction in line '{self._line}' must be positive")

		reac_parts = {
			"name": name.strip(),
			"rate": rate.strip(),
			"rate constant": float(k),
			"equation": eq.strip(),
			"coefficients": self._find_coefficients(eq.strip()),
			"num reactants": self._count_reactants(eq.strip()),
		}

		reaction.append(reac_parts)
		return reaction

	# reaction = irreversible, net, kinetics - arrhenius parameters A and Ea
	def _parse_irrev_net_AEa(self):
		# format of line: reaction name | rate law | arrhenius parameter A | arrhenius parameter Ea | reaction
		reaction = []

		parts = self._line.split("|")
		if len(parts) != 5:
			error = (
				"\n\nFor an irreversible net reaction controlled by arrhenius parameters A and Ea, the kinetic parameter(s) must be A and Ea (in J /mol)."
				"\ne.g. r1 | rate | 10e6 | 40000 | A -> B"
			)
			raise ValueError(self._bad_reaction_format(error))

		name, rate, A, Ea, eq = parts

		if float(A) <= 0 or float(Ea) <= 0:
			raise ValueError(f"Arrhenius parameters A and Ea in line '{self._line}' must both be positive")

		reac_parts = {
			"name": name.strip(),
			"rate": rate.strip(),
			"rate constant": self._calculate_k(float(A), float(Ea), self._temperature),
			"equation": eq.strip(),
			"coefficients": self._find_coefficients(eq.strip()),
			"A": A,
			"Ea": Ea,
			"num reactants": self._count_reactants(eq.strip()),
		}

		reaction.append(reac_parts)
		return reaction

	# reaction = reversible, elementary, kinetics = rate constant k
	def _parse_rev_elem_k(self):
		# format of line: reaction name | forward rate constant | backwards rate constant | reaction
		reaction = []

		parts = self._line.split("|")
		if len(parts) != 4:
			error = (
				"\n\nFor a reversible elementary reaction controlled by rate constants k, 2 kinetic parameters must be provided: rate constant kf of the forward reaction and rate constant kr of the reverse reaction"
				"\ne.g. r1 | 1e3 | 1e-2 | A <=> B"
			)
			raise ValueError(self._bad_reaction_format(error))

		name, kf, kr, eq = parts

		if float(kf) <= 0 or float(kr) <= 0:
			raise ValueError(f"Both rate constants for reaction in line '{self._line}' must be positive")

		lhs, rhs = eq.split("<=>")

		# forward reaction
		f_eq = f"{lhs} -> {rhs}"
		f_reac_parts = {
			"name": name.strip() + "_f",
			"rate constant": float(kf),
			"equation": f_eq.strip(),
			"coefficients": self._find_coefficients(f_eq.strip()),
			"num reactants": self._count_reactants(f_eq.strip()),
		}

		# reverse reaction
		r_eq = f"{rhs} -> {lhs}"
		r_reac_parts = {
			"name": name.strip() + "_r",
			"rate constant": float(kr),
			"equation": r_eq.strip(),
			"coefficients": self._find_coefficients(r_eq.strip()),
			"num reactants": self._count_reactants(r_eq.strip()),
		}

		reaction.append(f_reac_parts)
		reaction.append(r_reac_parts)
		return reaction 

	# reaction = reversible, elementary, kinetics = arrhenius parameters A and Ea
	def _parse_rev_elem_AEa(self):
		# format of line: reaction name | forwards arrhenius parameter A | forwards arrhenius parameter Ea | reverse arrhenius parameter A | reverse arrhenius parameter Ea | reaction
		reaction = []

		parts = self._line.split("|")
		if len(parts) != 6:
			error = (
				"\n\nFor a reversible elementary reaction controlled by arrhenius parameters A and Ea, the kinetic parameter(s) are two sets of arrhenius parameters A and Ea (in J /mol) must be provided, for the forward and reverse reactions."
				"\ne.g. r1 | 10e6 | 40000 | 10e9 | 80000 | A <=> B"
			)
			raise ValueError(self._bad_reaction_format(error))

		name, A_f, Ea_f, A_r, Ea_r, eq = parts

		if float(A_f) <= 0 or float(Ea_f) <= 0 or float(A_r) <= 0 or float(Ea_r) <= 0:
			raise ValueError(f"Arrhenius parameters A and Ea for the forward and reverse reactions in line '{self._line}' must all be positive")

		lhs, rhs = eq.split("<=>")

		# forward reaction
		f_eq = f"{lhs} -> {rhs}"
		f_reac_parts = {
			"name": name.strip() + "_f",
			"rate constant": self._calculate_k(float(A_f), float(Ea_f), self._temperature),
			"equation": f_eq.strip(),
			"coefficients": self._find_coefficients(f_eq.strip()),
			"A": A_f,
			"Ea": Ea_f,
			"num reactants": self._count_reactants(f_eq.strip()),
		}

		# reverse reaction
		r_eq = f"{rhs} -> {lhs}"
		r_reac_parts = {
			"name": name.strip() + "_r",
			"rate constant": self._calculate_k(float(A_r), float(Ea_r), self._temperature),
			"equation": r_eq.strip(),
			"coefficients": self._find_coefficients(r_eq.strip()),
			"A": A_r,
			"Ea": Ea_r,
			"num reactants": self._count_reactants(r_eq.strip()),
		}

		reaction.append(f_reac_parts)
		reaction.append(r_reac_parts)
		return reaction 

	# reaction = reversible, net, kinetics = rate constant k
	def _parse_rev_net_k(self):
		# format of line: reaction name | forward rate law | reverse rate law | forwards rate constant | reverse rate constant | reaction
		reaction = []

		parts = self._line.split("|")
		if len(parts) != 6:
			error = (
				"\n\nFor a reversible net reaction contorlled by rate constants k, make sure the reaction line is formatted as above."
				"\n2 rate laws must be provided, one for the forwards reaction and one for the reverse reaction. Rate laws must be of the format: k*[A]^2, where k is the rate constant and A (and additional species) are species concentrations."
				"\n2 kinetic parameters must be provided: rate constant kf of the forward reaction and rate constant kr of the reverse reaction"
				"\ne.g. r1 | k*[A]^2 | k*[B] | 1e3 | 1e-2 | A <=> B"
			)
			raise ValueError(self._bad_reaction_format(error))

		name, f_rate, r_rate, kf, kr, eq = parts

		if float(A_f) <= 0 or float(Ea_f) <= 0 or float(A_r) <= 0 or float(Ea_r) <= 0:
			raise ValueError(f"Arrhenius parameters A and Ea for the forward and reverse reactions in line '{self._line}' must all be positive")

		lhs, rhs = eq.split("<=>")

		# forward reaction
		f_eq = f"{lhs} -> {rhs}"
		f_reac_parts = {
			"name": name.strip() + "_f",
			"rate": f_rate.strip(),
			"rate constant": float(kf),
			"equation": f_eq.strip(),
			"coefficients": self._find_coefficients(f_eq.strip()),
			"num reactants": self._count_reactants(f_eq.strip()),
		}

		# reverse reaction
		r_eq = f"{rhs} -> {lhs}"
		r_reac_parts = {
			"name": name.strip() + "_r",
			"rate": r_rate.strip(),
			"rate constant": float(kr),
			"equation": r_eq.strip(),
			"coefficients": self._find_coefficients(r_eq.strip()),
			"num reactants": self._count_reactants(r_eq.strip()),
		}

		reaction.append(f_reac_parts)
		reaction.append(r_reac_parts)
		return reaction 

	# reaction = reversible, net, kinetics = arrhenius parameters A and Ea
	def _parse_rev_net_AEa(self):
		# format of line: reaction name | forward rate law | reverse rate law | forwards arrhenius parameter A | forwards arrhenius parameter Ea | reverse arrhenius parameter A | reverse arrhenius parameter Ea | reaction
		reaction = []

		parts = self._line.split("|")
		if len(parts) != 8:
			error = (
				"\n\nFor a reversible net reaction contorlled by arrhenius parameters A and Ea, make sure the reaction line is formatted as above."
				"\n2 rate laws must be provided, one for the forwards reaction and one for the reverse reaction. Rate laws must be of the format: k*[A]^2, where k is the rate constant and A (and additional species) are species concentrations."
				"\n4 kinetic parameters must be provided: arrhenius parameters A and Ea (in J /mol) of both the forward and reverse reactions, in the order: A (forward) | Ea (forward) | A (reverse) | Ea (reverse)."
				"\ne.g. r1 | k*[A]^2 | k*[B] | 10e6 | 40000 | 10e9 | 80000 | A <=> B"
			)
			raise ValueError(self._bad_reaction_format(error))

		name, f_rate, r_rate, A_f, Ea_f, A_r, Ea_r, eq = parts

		if float(A_f) <= 0 or float(Ea_f) <= 0 or float(A_r) <= 0 or float(Ea_r) <= 0:
			raise ValueError(f"Both rate constants for reaction in line '{self._line}' must be positive")

		lhs, rhs = eq.split("<=>")

		# forward reaction
		f_eq = f"{lhs} -> {rhs}"
		f_reac_parts = {
			"name": name.strip() + "_f",
			"rate": f_rate.strip(),
			"rate constant": self._calculate_k(float(A_f), float(Ea_f), self._temperature),
			"equation": f_eq.strip(),
			"coefficients": self._find_coefficients(f_eq.strip()),
			"A": A_f,
			"Ea": Ea_f,
			"num reactants": self._count_reactants(f_eq.strip()),
		}

		# reverse reaction
		r_eq = f"{rhs} -> {lhs}"
		r_reac_parts = {
			"name": name.strip() + "_r",
			"rate": r_rate.strip(),
			"rate constant": self._calculate_k(float(A_r), float(Ea_r), self._temperature),
			"equation": r_eq.strip(),
			"coefficients": self._find_coefficients(r_eq.strip()),
			"A": A_r,
			"Ea": Ea_r,
			"num reactants": self._count_reactants(r_eq.strip()),
		}

		reaction.append(f_reac_parts)
		reaction.append(r_reac_parts)
		return reaction 






