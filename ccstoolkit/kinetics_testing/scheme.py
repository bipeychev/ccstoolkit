
## SCHEME parses the input data, and produces the stoichiometric matrix, and the production rates required for simulation.

import numpy as np
from sympy import Matrix, lcm
from scipy import constants
import re
import ast
import os

from reaction_rate import elementary_reaction_rate, net_reaction_rate
from parse_reactions import PARSE_REACTION

class SCHEME():

	def __init__(self, inputData=None):
		self._inputData = inputData

		self._reactions = []
		self._species = []
		self._reactor_type = None
		self._qss_true = False
		self._reaction_time = None
		self._residence_time = None
		self._N_grid = None
		self._reaction_type = None
		self._temperature = None
		self._data_type = None
		self._parameter_type = None
		self._data = None

		if self._inputData is None:
			raise ValueError("Either an input file or input text must be supplied.")
		if os.path.isfile(inputData):
			self._load_file(inputData)
		else:
			self._load_text(inputData)

		if not self._reactions:
			raise ValueError("No reactions specified in input.")

		self._build_arrays()
		self._stoichiometric_matrix()

	def _load_file(self, filename):
		with open(filename, "r") as file:
			inputText = file.read()

		self._load_text(inputText)

	def _load_text(self, inputText):
		lines = [line.strip() for line in inputText.splitlines() if line.strip()]
		self._reaction_lines = []
		self._species_lines = []
		mode = None

		for line in lines:
			if line.startswith("#"):
				continue
			if line == "[CONDITIONS]":
				mode = "conditions"
				continue
			if line == "[REACTIONS]":
				mode = "reactions"
				continue
			if line == "[SPECIES]":
				mode = "species"
				continue
			if line == "[EXPERIMENTAL]":
				mode = "experimental"
				continue

			if mode == "reactions":
				self._reaction_lines.append(line)

			if mode == "species":
				self._species_lines.append(line)

			if mode == "conditions":
				if line.startswith("REACTOR TYPE"):
					self._reactor_type = (line.split("=")[1]).strip().lower()
					continue
				elif line.startswith("QSS"):
					self._qss_true = line.split("=")[1].strip().lower() == "true"
					continue
				elif line.startswith("REACTION TIME"):
					self._reaction_time = float(line.split("=")[1])
					continue
				elif line.startswith("RESIDENCE TIME"):
					self._residence_time = float(line.split("=")[1])
					continue
				elif line.startswith("GRID POINTS"):
					self._N_grid = int(line.split("=")[1])
					continue
				elif line.startswith("REACTION TYPE"):
					self._reaction_type = (line.split("=")[1]).strip().lower()
					continue
				elif line.startswith("TEMPERATURE (K)"):
					self._temperature = float(line.split("=")[1])
					continue
				elif line.startswith("PARAMETER TYPE"):
					self._parameter_type = (line.split("=")[1]).strip().lower()
					continue

			if mode == "experimental":
				if line.startswith("DATA TYPE"):
					self._data_type = (line.split("=")[1]).strip().lower()
					continue
				elif line.startswith("DATA"):
					data = line.split("=", 1)[1].strip()
					if data:
						self._data = ast.literal_eval(line.split("=", 1)[1].strip())
						continue

		# check necessary data is provided and correct
		self._check_input()

		# parses species
		for line in self._species_lines:
			self._parse_species(line)

		# parses raections
		for line in self._reaction_lines:
			reaction = PARSE_REACTION(line, self._reaction_type, self._parameter_type, self._species, self._temperature, self._calculate_k)
			self._reactions.extend(reaction._reaction)

	# function for parsing species
	def _parse_species(self, line):
		parts = line.split("|")
		if len(parts) != 3:
			raise ValueError(f"Bad species line format: '{line}'. \nSee the sample input file for instructions")
		name, c0, type = parts

		species = {
			"name": name.strip(),
			"initial concentration": float(c0),
			"type": type.strip(),
		}
		self._species.append(species)

	def _build_arrays(self):
		self._mask_intermediates = [i['type']=='reactive intermediate' for i in self._species]
		self._mask_buffered = [i['type']=='buffered' for i in self._species]
		self._mask_stable = [i['type']=='stable' for i in self._species]
		self._mask_dynamic = ~np.array(self._mask_buffered)  # for non-qss simulations

		self._N_intermediates = sum(self._mask_intermediates)

		self._c_initial = np.array([i['initial concentration'] for i in self._species])
		self._c_initial_stable = self._c_initial[self._mask_stable]
		self._c_initial_buffered = self._c_initial[self._mask_buffered]
		self._c_initial_dynamic = self._c_initial[self._mask_dynamic]

		self._labels = np.array([i['name'] for i in self._species])
		self._labels_stable = self._labels[self._mask_stable]
		self._labels_buffered = self._labels[self._mask_buffered]
		self._labels_intermediates = self._labels[self._mask_intermediates]
		self._labels_dynamic = self._labels[self._mask_dynamic]
		
		self._reaction_names = [reaction["name"] for reaction in self._reactions]

	def _stoichiometric_matrix(self):
		n_species = len(self._species)
		n_reactions = len(self._reactions)

		S = np.zeros((n_species, n_reactions))

		for j, reaction in enumerate(self._reactions):
			S[:, j] = reaction["coefficients"]

		self._S = S

		lines = ["\nSTOICHIOMETRIC MATRIX:"]
		lines.append(f"{'':<8}" + "  ".join(f"{s['name']:>6}" for s in self._species))
		
		for j, reaction in enumerate(self._reactions):
			row = self._S[:, j]
			row_str = "  ".join(f"{x:6.2f}" for x in row)
			lines.append(f"{reaction['name'] + ':':<8}{row_str}")

		return "\n".join(lines)

	# get the rates of all reactions at a given c (c = list of concentrations)
	def _get_all_rates(self,c):
		all_rates = []
		if self._reaction_type == 'elementary':
			for reac in self._reactions:
				rate = elementary_reaction_rate(c, reac['coefficients'], reac['rate constant'])
				all_rates.append(rate)

		elif self._reaction_type == 'net':
			for reac in self._reactions:
				rate = net_reaction_rate(c, reac['rate'], reac['rate constant'], self._species)
				all_rates.append(rate)

		return all_rates

	# get the production rates for all species at a given c
	def _get_all_prod_rate(self,c):
		rates = self._get_all_rates(c)
		prod_rates = self._S @ rates
		return prod_rates

	def _get_rate_constants(self):
		rate_constants = np.array([reaction["rate constant"] for reaction in self._reactions])
		return rate_constants

	def _get_arrhenius_params(self):
		params = np.empty(2*len(self._reactions), dtype=float)
		A = np.array([reaction["A"] for reaction in self._reactions])
		Ea = np.array([reaction["Ea"] for reaction in self._reactions])
		params[::2] = A
		params[1::2] = Ea 
		return params

	def _set_rate_constants(self, ks):
		if len(ks) != len(self._reactions):
			raise ValueError(f"Number of updated rate constants ({len(ks)}) does not match number of original reactions ({len(self._reactions)}).")
		
		if np.any(np.array(ks) <= 0):
			raise ValueError(f"All updated rate constants must be positive")

		for reaction, ki in zip(self._reactions, ks):
			reaction["rate constant"] = ki

	def _check_input(self): # add other stuff here 
		# checking reactor type
		reactor_types = ['batch', 'pfr', 'cstr']
		if self._reactor_type == None:
			raise ValueError(f"Reactor type must be specified in conditions. Available reactor types are: {", ".join(reactor_types)}.")
		if self._reactor_type not in reactor_types:
			raise ValueError(f"Unknown reactor type '{self._reactor_type}'. Available reactor types are: {", ".join(reactor_types)}.")

		# checking necessary timing information is provided
		if self._reaction_time == None and self._reactor_type == "batch" or self._reaction_time == None and self._reactor_type == "cstr":
			raise ValueError(f"Reactor type '{self._reactor_type}' requires reaction time to be specified in input.")
		if self._residence_time == None and self._reactor_type == "pfr" or self._residence_time == None and self._reactor_type == "cstr":
			raise ValueError(f"Reactor type '{self._reactor_type}' requires residence time to be specified in input.")

		# checking type of reactions
		reaction_types = ['elementary', 'net']
		if self._reaction_type == None:
			raise ValueError(f"Reaction type must be specified in conditions. Available reaction types are: {", ".join(reaction_types)}.")
		if self._reaction_type not in reaction_types:
			raise ValueError(f"Unknown reaction type '{self._reaction_type}'. Available reaction types are: {", ".join(reaction_types)}.")

		# checking parameter type for optimising 
		parameter_types = ['rate constant', 'arrhenius']
		if self._parameter_type not in parameter_types:
			raise ValueError(f"Unknown parameter type: {self._parameter_type}. Available parameter types are: {", ".join(parameter_types)}. Please input this in [CONDITIONS].")

		# checking experimental data type (if data is provided)
		if self._data is not None:
			data_types = ['time course', 'final concentration']
			if self._data_type not in data_types:
				raise ValueError(f"Unknown data type: {self._data_type}. Available experimental data types are: {", ".join(data_types)}.")

		# checking temperature is provided (if parameter type = arrhenius)
		if self._parameter_type == "arrhenius":
			if self._temperature == None:
				raise ValueError("Please provide a simulation temperature. This is required for PARAMETER TYPE = arrhenius.")
			elif self._temperature <= 0:
				raise ValueError("Temperature must be provided in Kelvin and must be greater than 0 K.")

	# calculates rate constant k from A, Ea, T
	def _calculate_k(self, A, Ea, temperature):

		R = constants.gas_constant # J/mol/K
		T = temperature # K

		k = A * np.exp(-Ea/R/T)

		return k

	# takes an array of ln(A) and Ea, and T, and returns an array of ks
	def _calculate_ks(self, kinetic_parameters, temperature):
		# kinetic_parameters = [ln(A1), Ea1, ln(A2), Ea2, ...]

		rate_constants = []
		for i in range(0, len(kinetic_parameters), 2):
			A = np.exp(kinetic_parameters[i])
			Ea = kinetic_parameters[i+1]

			k = self._calculate_k(A, Ea, temperature)

			rate_constants.append(k)

		return np.asarray(rate_constants)