
## contains functions for fitting rate constants to experimental data

import numpy as np
import time
from scipy.optimize import least_squares
from fitting import FITTING

class RC_FITTING(FITTING):

	def __init__(self, inputData=None, optimised_reactions=None, initial_guesses=None):

		# gets required data from general FITTING class
		super().__init__(inputData=inputData, optimised_reactions=optimised_reactions, initial_guesses=initial_guesses)

	def _default_initial_guesses(self):
		# if initial guesses are not provided: uses the supplied rate constants in input
		rate_constants = self.scheme._get_rate_constants()
		initial_guesses = np.array([rate_constants[index] for index in self.optimised_reaction_index], dtype=float)

		return initial_guesses

	def _check_initial_guesses(self, initial_guesses):
		# checks given initial guesses are in the correct format 
		if len(init_initial_guesses) != len(self.optimised_reactions):
			raise ValueError(f"Number of initial guesses given for rate constants ({len(init_initial_guesses)}) must match number of optimised reactions ({len(self.optimised_reactions)}).")

		if np.any(init_initial_guesses <= 0):
			raise ValueError(f"Initial guesses given for rate constants must be positive")

	# residuals using rte constants k 
	def _rc_residuals(self, lg_rate_constants): 
		# lg_rate_constants = [ln(k1), ln(k2), ...]
		# conversion to normal space (rate_constants = [k1, k2, ...])
		rate_constants = np.exp(lg_rate_constants)

		# replace rate constants being adjusted
		new_rate_constants = self.scheme._get_rate_constants().copy()
		for index, k in zip(self.optimised_reaction_index, rate_constants):
			new_rate_constants[index] = k

		# calculate residuals for each experiment
		residuals = []
		for experiment in self.experimental_data:
			experiment_residuals = self._experiment_residuals(experiment, new_rate_constants)
			residuals.extend(experiment_residuals)

		return np.asarray(residuals)

	# fits rate constants k from experimental data 
	def fit_rate_constants(self, scan=False, method="trf", diff_step=0.01, ftol=1e-12, xtol=1e-12, gtol=1e-12, **kwargs):
		print(f"\nOptimising rate constants for reactions: {', '.join(self.optimised_reactions)}")

		if scan:
			# finds a reasonable initial guess for the system
			scan_guesses = self._k_log_scan()
			lg_guesses = np.log(scan_guesses) 
		else:
			# initial guesses in log space for optimisation
			lg_guesses = np.log(self.initial_guesses)

		# bounds
		lower_bound = np.log(np.full(len(self.optimised_reactions), 1e-10))
		upper_bound = np.log(np.full(len(self.optimised_reactions), 1e10))

		# minimise residuals 
		result = least_squares(self._rc_residuals, x0=lg_guesses, bounds=(lower_bound, upper_bound), args=(), method=method, diff_step=diff_step, ftol=ftol, xtol=xtol, gtol=gtol, **kwargs)

		# check result is well positioned
		if result.status == 1:
			raise ValueError(f"\nThe initial guess for the rate constant of reaction '{', '.join(self.optimised_reactions)}' is not optimal for this model. Either provide a different initial guess or pass scan=True to generate a more suitable starting point.")

		return result

	def optimise_rate_constants(self, scan=False, method="trf", diff_step=0.01, ftol=1e-12, xtol=1e-12, gtol=1e-12, **kwargs):
		optimise_start_time = time.perf_counter()

		# optimises rate constants based on experimental data 
		result = self.fit_rate_constants(method=method, diff_step=diff_step, scan=scan, ftol=ftol, xtol=xtol, gtol=gtol, **kwargs)

		fitted_rate_constants = np.exp(result.x)
		dk_lower, dk_upper = self._rate_constant_uncertainty(result)

		optimise_end_time = time.perf_counter()
		optimisation_time = optimise_end_time - optimise_start_time

		#store optimised rate constants + corresponding uncertainty
		result.fitted_rate_constants = fitted_rate_constants
		result.uncertainties = dk_lower, dk_upper
		result.optimisation_time = optimisation_time

		# print results
		print(f"\nOptimisation time: {optimisation_time:.4} s")
		print(f"Optimised result SSR: {np.sum(result.fun**2):.4}")
		print(f"\nOptimisation successful. The optimised rate constants (with upper and lower uncertainties) are:")

		for i, reac in enumerate(self.optimised_reactions):
			units = self._determine_units(reac)
			print(f"k ({reac}): {fitted_rate_constants[i]:.5g} - {dk_lower[i]:.5g} + {dk_upper[i]:.5g} {units}")

		return result

	# scans log space from 1e-10 to 1e10 for a good starting guess for the rate constants - one at a time, not multivariable
	def _k_log_scan(self):
		print(f"\nPerforming scan for reasonable starting point.")

		scan_time_start = time.perf_counter()

		# create log spaced out k values
		scan_ks = np.logspace(-10, 15, 26)

		# scan each value to find the lowest residual
		best_guesses = self.initial_guesses.copy()
		for i, reac in enumerate(self.optimised_reactions):
			# test each value of k
			scan_ssr = []
			for k in scan_ks:
				test_guesses = best_guesses.copy()

				# replace k for this reaction
				test_guesses[i] = k

				# calculate residuals
				residuals = self._rc_residuals(np.log(test_guesses))
				scan_ssr.append(np.sum(residuals**2))

			# find the k with the lowest SSR
			best_k = scan_ks[np.argmin(scan_ssr)]
			best_guesses[i] = best_k

		scan_time_end = time.perf_counter()
		scan_time = scan_time_end - scan_time_start

		print(f"\nReasonable starting points for the rate constants of these reactions:")
		for reac, k in zip(self.optimised_reactions, best_guesses):
			units = self._determine_units(reac)
			print(f"{reac}: {k:.5g} {units}")
		print(f"Scan time: {scan_time:.4} s")

		return best_guesses

	def _rate_constant_uncertainty(self, result):
		# jacobian in log(k) space 
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

		# log uncertainties
		sigma_k = np.sqrt(np.diag(covariance)) # in log space

		# conversion to normal space 
		fitted_k = np.exp(result.x) # [k1, k2, ...]
		dk_lower = np.abs(fitted_k-np.exp(result.x-sigma_k)) # lower uncertainty
		dk_upper = np.abs(fitted_k-np.exp(result.x+sigma_k)) # upper uncertainty

		return dk_lower, dk_upper

