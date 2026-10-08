
## class to fit kinetic parameters (either k or A & Ea) from experimental data

import numpy as np
import matplotlib.pyplot as plt
from simulate import simulate_scheme
from scheme import SCHEME

class FITTING():

	def __init__(self, inputData=None, optimised_reactions=None, initial_guesses=None):
		if inputData is None:
			raise ValueError(f"Input data must be provided")
		if isinstance(inputData, str):
			inputData = [inputData]

		self.inputData = inputData
		self.schemes = [SCHEME(input) for input in inputData]
		self.scheme = self.schemes[0]

		# reaction + species info
		self.all_reactions = self.scheme._reaction_names
		self.all_species = [sp['name'] for sp in self.scheme._species]

		# dictionary of experimental data
		self.experimental_data = self._parse_experimental_data()

		# list of all species with experimental data
		self.expt_species = sorted({sp for experiment in self.experimental_data for sp in experiment['data']})

		# list of all temperatures
		self._all_temperatures = [experiment["temperature"] for experiment in self.experimental_data]

		# list of reactions with kinetic parameters to optimise
		self.optimised_reactions = self._set_optimised_reactions(optimised_reactions)
		self.optimised_reaction_index = [self.all_reactions.index(reaction) for reaction in self.optimised_reactions]

		# array of initial guesses
		self.initial_guesses = self._set_initial_guesses(initial_guesses)

	def _parse_experimental_data(self):

		# check all experiments use the same reaction scheme
		for i, scheme in enumerate(self.schemes[1:], start=2):
			reactions = scheme._reaction_names
			species = [sp['name'] for sp in scheme._species]

			if reactions != self.all_reactions:
				raise ValueError(f"Experiment {i} has a different reaction scheme to experiment 1.")
			if species != self.all_species:
				raise ValueError(f"Experiment {i} has different species to experiment 1.")

		# create dictionary of expeirmental data
		experimental_data = []
		for i, scheme in enumerate(self.schemes, start=1):
			if scheme._data_type is None:
				raise ValueError(f"Experiment {i} is missing DATA TYPE. This is required for kinetic parameter optimisation. Options are 'time course' (series of time and concetration values provided) or 'final conc' (one time and concentration value).")
			if scheme._data is None:
				raise ValueError(f"Experiment {i} is missing DATA. This is required for kinetic parameter optimisation. Provide data relevant to the provided DATA TYPE for 1 or more species.")

			experiment = {
				"name": f"Experiment {i}",
				"input": scheme._inputData,
				"type": scheme._data_type.replace(" ", "_"),
				"data": scheme._data,
				"temperature": scheme._temperature,
			}
			experimental_data.append(experiment)

		return experimental_data

	def _set_optimised_reactions(self, init_optimised_reactions):
		# if optimised reactions are not provided: uses check_rc_dependence to ascertain reaction rate constants dependent on the experimental data 
		if init_optimised_reactions is None:
			optimised_reactions = sorted(self.check_rc_dependence())

		else:
			# check inputted reactions are in the correct format
			if not isinstance(init_optimised_reactions, list):
				raise TypeError("Supplied optimised reactions must be a list of strings, e.g. ['r1', 'r2'].")
			if not all(isinstance(reaction, str) for reaction in init_optimised_reactions):
				raise TypeError("Supplied optimised reactions must be a list of strings, e.g. ['r1', 'r2'].")

			# check that requested reactions exist in the scheme
			for reac in init_optimised_reactions:
				if reac not in self.all_reactions:
					raise ValueError(f"Reaction specified for optimisation '{reac}' not found. Available reactions: {self.all_reactions}")

			optimised_reactions = init_optimised_reactions

		return optimised_reactions

	def _set_initial_guesses(self, init_initial_guesses):
		if init_initial_guesses is None:
			initial_guesses = self._default_initial_guesses()

		else:
			self._check_initial_guesses(init_initial_guesses)
			initial_guesses = np.asarray(init_initial_guesses, dtype=float)

		return initial_guesses

	# calulates the residual for one set of experimental data
	def _experiment_residuals(self, experiment, new_rate_constants):

		# simulates scheme with new rate constants
		simulated_system, _ = simulate_scheme(experiment['input'], rate_constants=new_rate_constants)

		# calculate residuals for each experiment
		residuals = []
		for sp, data in experiment['data'].items():

			# time course data
			if experiment['type'] == "time_course":
				# simulation data
				simulation_time = simulated_system.result[sp]['x']
				simulation_values = simulated_system.result[sp]['y']

				# experiment data
				experimental_time = np.asarray(data['time'], dtype=float)
				experimental_values = np.asarray(data['values'], dtype=float)
				if len(experimental_time) != len(experimental_values):
					raise ValueError(f"Experimental data {experiment['name']} has {len(experimental_time)} time points but {len(experimental_values)} experimental data points for species {sp}.")

				# interpolate simulated data to match experimental time points
				simulated_values = np.interp(experimental_time, simulation_time, simulation_values)

				species_residuals = (simulated_values - experimental_values)
				residuals.extend(species_residuals)

			# final concentration data
			elif experiment['type'] == "final_conc":
				# simulated value & experimental data point
				simulated_value = simulated_system.result[sp]['y'][-1]
				experimental_value = float(data)

				residual = simulated_value - experimental_value
				residuals.append(residual)

		return np.asarray(residuals)

	# determines which reaction rate constants affect the experimental data. doesnt prove that A and Ea are seperately identifiable
	def check_rc_dependence(self):
		# get provided rate constants (k). for given arrhenius parameters k is calculated based on A, Ea, T
		original_rate_constants = self.scheme._get_rate_constants()

		print(f"\nTesting the rate constant dependence across all experiments. A rate constant is considered arbitrary if a significant change does not affect the concentration.")

		dependent_reactions = set()

		# tests each reaction rate constant
		for index, reaction_name in enumerate(self.all_reactions):
			dependence = False
			original_k = original_rate_constants[index]

			# change rate constant by 1e-3, 1e3
			test_values = [original_k / 1000, original_k * 1000]
			for test_k in test_values:
				test_rate_constants = original_rate_constants.copy()
				test_rate_constants[index] = test_k

				# test each experiment
				for experiment in self.experimental_data:
					# baseline simulation
					baseline_system, _ = simulate_scheme(experiment['input'], rate_constants=original_rate_constants)

					# perturbed simulation
					test_system, _ = simulate_scheme(experiment['input'], rate_constants=test_rate_constants)

					# check the species with experimental data 
					for sp, experimental_data in experiment['data'].items():
						if experiment['type'] == "time_course":
							experimental_time = np.asarray(experimental_data['time'], dtype=float)

							# interpolate baseline + test simulations to experimental times
							baseline_data = np.interp(experimental_time, baseline_system.result[sp]['x'], baseline_system.result[sp]['y'])
							test_data = np.interp(experimental_time, test_system.result[sp]['x'], test_system.result[sp]['y'])

						elif experiment['type'] == "final_conc":
							baseline_data = baseline_system.result[sp]['y'][-1]
							test_data = test_system.result[sp]['y'][-1]

						percentage_change = np.abs(100 * (baseline_data - test_data) / baseline_data)

						if np.any(percentage_change >= 1):
							dependence = True
							break 
					if dependence:
						break
				if dependence:
					break
			if dependence:
				dependent_reactions.add(reaction_name)

		print(f"\nReactions with non-arbitrary rate constants: {', '.join(dependent_reactions)}")
		return dependent_reactions

	# determines units for k or A
	def _determine_units(self, reac):
		for reaction in self.scheme._reactions:
			if reaction["name"] == reac:
				if reaction["num reactants"] == 1:
					units = "/s"
				elif reaction["num reactants"] == 2:
					units = "/mM/s"
				else:
					units = "units in the equivalent of /mM/s dependent on number of reagents"

		return units


