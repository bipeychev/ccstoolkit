
## contains functions for fitting arrhenius parameters to experimental data

import numpy as np
import time
from scipy.optimize import least_squares, differential_evolution
from fitting import FITTING

class ARRHENIUS_FITTING(FITTING):

	def __init__(self, inputData=None, optimised_reactions=None, initial_guesses=None):

		# gets required data from general FITTING class
		super().__init__(inputData=inputData, optimised_reactions=optimised_reactions, initial_guesses=initial_guesses)

	def _default_initial_guesses(self):
		# if initial guesses are not provided: uses the supplied rate constants in input
		arrhenius_params = self.scheme._get_arrhenius_params()
		initial_guesses = np.concatenate([arrhenius_params[2*index:2*index+2] for index in self.optimised_reaction_index], dtype=float)

		return initial_guesses

	def _check_initial_guesses(self, init_initial_guesses):
		# checks given initial guesses are in the correct format [ln(A1), Ea1, ln(A2), Ea2,...]
		if len(init_initial_guesses) != 2*len(self.optimised_reactions):
			raise ValueError(f"Expected 2 initial guesses (A and Ea) for each optimised reaction (total {2*len(self.optimised_reactions)}), but recieved {len(init_initial_guesses)}.")

		if np.any(init_initial_guesses[::2] <= 0):
			raise ValueError("Initial guesses for A must be > 0")
		if np.any(init_initial_guesses[1::2] < 0):
			raise ValueError("Initial guesses for Ea must be >= 0")

	# residuals for optimising arrhenius parameters A and Ea
	def _arrhenius_residuals(self, lg_arrhenius_vals):
		# lg_arrhenius_vals: [ln(A1), Ea1, ln(A2), Ea2, ...]

		# calculates residuals for each experiment
		residuals = []
		for experiment in self.experimental_data:

			# replace rate constants being adjusted - temperature dependent
			rate_constants = self.scheme._get_rate_constants().copy()
			temperature = experiment['temperature']
			new_rate_constants = self.scheme._calculate_ks(lg_arrhenius_vals, temperature)

			for index, k in zip(self.optimised_reaction_index, new_rate_constants):
				rate_constants[index] = k

			# calculate residual
			experiment_residuals = self._experiment_residuals(experiment, rate_constants)
			residuals.extend(experiment_residuals)

		return np.asarray(residuals)

	# fits arrhenius A and Ea values
	def fit_arrhenius(self, scan=False, global_search=False, method="trf", diff_step=None, x_scale='jac', ftol=1e-12, xtol=None, gtol=1e-12, **kwargs):
		print(f"\nOptimising arrhenius parameters A and Ea for reactions: {', '.join(self.optimised_reactions)}")

		# bounds 
		lower_bound = np.empty(2*len(self.optimised_reactions))
		upper_bound = np.empty(2*len(self.optimised_reactions))

		A_min = 1e-10
		A_max = 1e30
		lower_bound[::2] = np.log(A_min)
		upper_bound[::2] = np.log(A_max)

		Ea_min = 0
		Ea_max = 1e6
		lower_bound[1::2] = Ea_min
		upper_bound[1::2] = Ea_max

		if global_search:
			lgA_guesses = self._global_search(lower_bound, upper_bound)
		elif scan:
			# finds a reasonable initial guess for the system
			scan_guesses = self._arrhenius_log_scan()
			lgA_guesses = self._log_convert(scan_guesses)
		else:
			# conversion from [A, Ea] to [ln(A), Ea]
			lgA_guesses = self._log_convert(self.initial_guesses)

		# minimise residuals 
		result = least_squares(self._arrhenius_residuals, x0=lgA_guesses, bounds=(lower_bound, upper_bound), method=method, diff_step=diff_step, ftol=ftol, xtol=xtol, gtol=gtol, **kwargs)

		### maybe change this - look at other func
		### add global optimisation aswell as local 
		if not result.success:
			raise ValueError(f"Optimisation did not converge. Solver message: {result.message}")

		return result

	def optimise_arrhenius(self, scan=False, global_search=False, method="trf", diff_step=0.01, ftol=1e-12, xtol=1e-12, gtol=1e-12, **kwargs):
		optimise_start_time = time.perf_counter() 

		# optimises A and Ea based on experimental data 
		result = self.fit_arrhenius(method=method, diff_step=diff_step, scan=scan, global_search=global_search, ftol=ftol, xtol=xtol, gtol=gtol, **kwargs)

		fitted_parameters = result.x # format [ln(A1), Ea1, ln(A2), Ea2, ...]
		lower, upper = self._arrhenius_uncertainty(result)

		optimise_end_time = time.perf_counter()
		optimisation_time = optimise_end_time - optimise_start_time

		# store results and corresponding uncertainties
		result.fitted_parameters = fitted_parameters
		result.fitted_A = np.exp(fitted_parameters[::2]) # [A1, A2,..]
		result.fitted_Ea = fitted_parameters[1::2] 
		result.uncertainties = lower, upper
		result.optimisation_time = optimisation_time

		# print results
		print(f"\nOptimisation time: {optimisation_time:.4} s")
		print(f"Optimised SSR: {np.sum(result.fun**2):.4}")
		print(f"\nOptimisation successful. The optimised arrhenius parameters (with upper and lower uncertainties) are:")
		for i, reac in enumerate(self.optimised_reactions):
			units = self._determine_units(reac)
			print(f"A ({reac}): {result.fitted_A[i]:.5g} - {lower[2*i]:.5g} + {upper[2*i]:.5g} {units}") 
			print(f"Ea ({reac}): {result.fitted_Ea[i]:.5g} - {lower[2*i+1]:.5g} + {upper[2*i+1]:.5g} J /mol")

		return result

	# scans log space for best initial guess of A and Ea
	def _arrhenius_log_scan(self):
		print(f"\nPerforming scan for reasonable starting point.")

		scan_time_start = time.perf_counter()

		# scans for a good starting guess for A and Ea
		# total scans = 15^2 = 225
		scan_As = np.logspace(-10, 30, 21) # scan A logarithmically
		scan_Eas = np.linspace(0, 1e6, 21) # scan Ea linearly

		# scan each value to find the lowest residual
		best_guesses = self.initial_guesses.copy() # [A, Ea]
		for i, reac in enumerate(self.optimised_reactions):
			# test each value of A and Ea 
			scan_ssr = []
			scan_guesses = []

			for A in scan_As:
				for Ea in scan_Eas:
					# test validity of resulting k at chosen A and Ea values for experimental temperatures
					valid = True
					for temp in self._all_temperatures:
						k = self.scheme._calculate_k(A, Ea, temp)
						if k <= 1e-10 or k >= 1e15:
							valid = False
							break
					if not valid:
						continue

					test_guesses = best_guesses.copy()

					# replace A and Ea for this reaction
					test_guesses[2*i] = A
					test_guesses[2*i+1] = Ea

					# conversion from [A, Ea] to [ln(A), Ea]
					lg_test_guesses = self._log_convert(test_guesses)

					# calculate residuals
					residuals = self._arrhenius_residuals(lg_test_guesses)
					scan_ssr.append(np.sum(residuals**2))
					scan_guesses.append((A,Ea))

			if len(scan_ssr) == 0:
				raise ValueError(f"No valid A/Ea combinations were found for {reac}. Check the A and Ea scan ranges.")

			# find A and Ea combination with the lowest SSR
			best_A, best_Ea = scan_guesses[np.argmin(scan_ssr)]
			best_guesses[2*i] = best_A
			best_guesses[2*i+1] = best_Ea

		scan_time_end = time.perf_counter()
		scan_time = scan_time_end - scan_time_start

		print(f"\nReasonable starting points for the arrhenius parameters of these reactions:")
		for reac, A, Ea in zip(self.optimised_reactions, best_guesses[::2], best_guesses[1::2]):
			units = self._determine_units(reac)
			print(f"{reac}: A = {A:.5g} {units}, Ea = {Ea:.5g} J /mol")
		print(f"Scan time: {scan_time:.4} s")

		return best_guesses

	def _arrhenius_uncertainty(self, result):
		# optimiser result is in [ln(A), Ea] space
		J = result.jac

		# number of observations & fitted parameters
		n_data = len(result.fun)
		n_parameters = len(result.x)
		degrees_of_freedom = n_data - n_parameters

		#optimised deviation
		if degrees_of_freedom > 0:
			opt_dev = np.sum((result.fun)**2) / degrees_of_freedom
		else:
			opt_dev = np.sum((result.fun)**2) / n_data

		# covariance matrix, K = opt_dev * 1/2 * H^-1
		# for least squares: H = 2 J^T * J 
		# K = (J^T * J)^-1
		JTJ = J.T @ J

		covariance = opt_dev * np.linalg.inv(JTJ)

		# uncertainties
		sigma = np.sqrt(np.diag(covariance))

		fitted_params = result.x # [ln(A1), Ea1, ln(A2), Ea2, ...]
		fitted_A = np.exp(fitted_params[::2]) # [A1, A2,..]
		fitted_Ea = fitted_params[1::2]  # [Ea1, Ea2, ...]

		sigma_lnA = sigma[::2]
		A_lower = np.abs(fitted_A - np.exp(fitted_params[::2] - sigma_lnA))
		A_upper = np.abs(fitted_A - np.exp(fitted_params[::2] + sigma_lnA))

		sigma_Ea = sigma[1::2]
		Ea_lower = sigma_Ea.copy()
		Ea_upper = sigma_Ea.copy()

		# store all uncertainties
		lower = np.empty_like(result.x)
		upper = np.empty_like(result.x)
		lower[::2] = A_lower
		lower[1::2] = Ea_lower
		upper[::2] = A_upper
		upper[1::2] = Ea_upper

		return lower, upper 

	def _log_convert(self, params):
		# params: [A1, Ea1, A2, Ea2, ...]
		# lg_params: [ln(A1), Ea1, ln(A2), Ea2, ...]

		# conversion from [A, Ea] to [ln(A), Ea]
		lg_params = params.copy()
		lg_params[::2] = np.log(params[::2])

		return lg_params

	def _global_search(self, lower_bound, upper_bound, seed=None, maxiter=50, popsize=10):
		print("\nPerforming global search (differential evolution)")

		def ssr(params):
			# params = [ln(A1), Ea1, ln(A2), Ea2, ...]
			try:
				res = self._arrhenius_residuals(params)
			except Exception:
				return 1e30 # solver failure - reject this point
			if not np.all(np.isfinite(res)):
				return 1e30
			return np.sum(res**2)

		result = differential_evolution(ssr, bounds=list(zip(lower_bound, upper_bound)), seed=seed, maxiter=maxiter, popsize=popsize, init="sobol", tol=1e-6, polish=False,)

		print(f"Global search SSR: {result.fun:.4g}")
		return result.x

# leave maxiter and popsize as default.
# test difference with least squares after/ not