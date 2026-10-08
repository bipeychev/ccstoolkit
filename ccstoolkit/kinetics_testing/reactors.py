
## contains the classes for the non steady-state reactor simulation codes

import numpy as np
from scipy.integrate import solve_ivp

class BATCH_REACTOR():

	def __init__(self, scheme):
		self.scheme = scheme

	def _sode(self, t, c_dynamic): # swapped round for diff integrator
		# list of all concentrations
		c = np.zeros(len(self.scheme._species))
		c[self.scheme._mask_buffered] = self.scheme._c_initial_buffered
		c[self.scheme._mask_dynamic] = c_dynamic

		# list of production rates for each species
		production_rates = self.scheme._get_all_prod_rate(c)
		
		# return production rates for all non-buffered species
		return production_rates[self.scheme._mask_dynamic]

	def simulate(self, t_max, N_grid=None):
		if N_grid == None:
			raise ValueError("Please provide a value for N_grid. Suggested: 1000")

		t = np.linspace(0, t_max, N_grid)
		sol = solve_ivp(self._sode, (0, t_max), self.scheme._c_initial_dynamic, method="BDF", t_eval=t, rtol=1e-8, atol=1e-12)

		c_buffered = {} #do i need this???
		for index, species in enumerate(self.scheme._labels_buffered):
			c_buffered[species] = {
				"x": t,
				"y": np.array([self.scheme._c_initial_buffered[index]] * N_grid)
			}

		c_dynamic = {}
		for species, conc in zip(self.scheme._labels_dynamic, sol.y):
			c_dynamic[species] = {
				"x": t,
				"y": conc
			}

		self.result = c_buffered | c_dynamic
		#print("Batch reactor simulation finished successfully")

class PLUG_FLOW_REACTOR():
	# differs from batch only by use of tau (residence time)

	def __init__(self, scheme):
		self.scheme = scheme

	def _sode(self, tau, c_dynamic): 
		c = np.zeros(len(self.scheme._species))
		c[self.scheme._mask_buffered] = self.scheme._c_initial_buffered
		c[self.scheme._mask_dynamic] = c_dynamic

		# list of production rates for each species
		production_rates = self.scheme._get_all_prod_rate(c)
		
		# return production rates for all non-buffered species
		return production_rates[self.scheme._mask_dynamic]

	def simulate(self, tau_max, N_grid=None):
		if N_grid == None:
			raise ValueError("Please provide a value for N_grid. Suggested: 1000")

		tau = np.linspace(0, tau_max, N_grid)
		#sol = odeint(self._sode, self.scheme._c_initial_dynamic, t)
		sol = solve_ivp(self._sode, (0, tau_max), self.scheme._c_initial_dynamic, method="BDF", t_eval=tau, rtol=1e-8, atol=1e-12)

		c_buffered = {} # not sure if necessary
		for index, species in enumerate(self.scheme._labels_buffered):
			c_buffered[species] = {
				"x": tau,
				"y": np.array([self.scheme._c_initial_buffered[index]] * N_grid)
			}

		c_dynamic = {}
		for species, conc in zip(self.scheme._labels_dynamic, sol.y):
			c_dynamic[species] = {
				"x": tau,
				"y": conc
			}

		self.result = c_buffered | c_dynamic
		#print("PFR simulation finished successfully")

class CSTR_REACTOR():
	def __init__(self, scheme):
		self.scheme = scheme

	def _sode(self, t, c_dynamic, tau): 
		c = np.zeros(len(self.scheme._species))
		c[self.scheme._mask_buffered] = self.scheme._c_initial_buffered
		c[self.scheme._mask_dynamic] = c_dynamic

		# list of production rates for each species
		production_rates = self.scheme._get_all_prod_rate(c)

		# differs here - cstr rate = production rate + (C_0 - C)/tau
		flow_term = np.zeros(len(self.scheme._species))
		flow_term[self.scheme._mask_dynamic] =  (self.scheme._c_initial_dynamic - c_dynamic)/tau

		# return production rates for all non-buffered species
		cstr_rates = production_rates[self.scheme._mask_dynamic] + flow_term[self.scheme._mask_dynamic]
		return cstr_rates

	def simulate(self, t_max, tau, N_grid=None):
		if N_grid == None:
			raise ValueError("Please provide a value for N_grid. Suggested: 1000")

		t = np.linspace(0, t_max, N_grid)
		sol = solve_ivp(self._sode, (0, t_max), self.scheme._c_initial_dynamic, method="BDF", t_eval=t, args=(tau,), rtol=1e-8, atol=1e-12)

		c_buffered = {} #do i need this???
		for index, species in enumerate(self.scheme._labels_buffered):
			c_buffered[species] = {
				"x": t,
				"y": np.array([self.scheme._c_initial_buffered[index]] * N_grid)
			}

		c_dynamic = {}
		for species, conc in zip(self.scheme._labels_dynamic, sol.y):
			c_dynamic[species] = {
				"x": t,
				"y": conc
			}

		self.result = c_buffered | c_dynamic
		#print("CSTR simulation finished sucessfully")
