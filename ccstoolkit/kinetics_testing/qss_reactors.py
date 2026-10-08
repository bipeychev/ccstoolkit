
## contains the classes for the steady-state reactor simulation codes

import numpy as np
from scipy.optimize import root, minimize
from scipy.integrate import odeint

class BATCH_REACTOR_QSS():

	def __init__(self, scheme):
		self.scheme = scheme

	#System of ordinary equations for the reactive intermediates at steady-state
	def _soe(self,c_intermediates,c_stable):
		
		#List of all concentrations
		c = np.zeros(len(self.scheme._species))
		c[self.scheme._mask_buffered] = self.scheme._c_initial_buffered
		c[self.scheme._mask_intermediates] = c_intermediates
		c[self.scheme._mask_stable] = c_stable

		#List of production rates for each species
		production_rates = self.scheme._get_all_prod_rate(c)
		
		#Return only the production rates for the reactive intermediates
		return production_rates[self.scheme._mask_intermediates]
	
	#Get the steady-state concentration of the intermediates given the concentrations of the other species
	def _get_intermediates(self,c_stable):
		#Initial guess for the concentrations of the reactive intermediates
		lgc_int_0 = np.array([-6]*self.scheme._N_intermediates)

		#Wrapper for the functions to be minimized
		min_function = lambda lgc_int: np.sum([i**2 for i in self._soe(10**lgc_int,c_stable)])
		
		#Find the root of _soe, i.e. the concentrations of the reactive intermediates that minimize the sum of squares
		sol = minimize(min_function, lgc_int_0, method='Nelder-Mead', tol=1e-8, options={'maxiter': 1000, 'xatol': 1e-8, 'fatol': 1e-8, 'adaptive': True})
		#'adaptive': True is very important!

		if not sol.success:
			print("SOE for the reactive intermediates did not solve.")
			print(f"concentrations of stable species: {c_stable}")
			print(sol)
			return [np.Nan]*self.scheme._N_intermediates

		return 10**sol.x
	
	#System of ordinary differential equations for the non-buffered, non-steady state intermediates
	def _sode(self,c_stable, t):																	  #~40ms

		#Solve for the reactive intermediates
		c_intermediates = self._get_intermediates(c_stable)
		
		#List of all concentrations
		c = np.zeros(len(self.scheme._species))
		c[self.scheme._mask_buffered] = self.scheme._c_initial_buffered
		c[self.scheme._mask_intermediates] = c_intermediates
		c[self.scheme._mask_stable] = c_stable

		#List of production rates for each species
		production_rates = self.scheme._get_all_prod_rate(c)
		
		#print(production_rates)
		#Return only the production rates for the non-buffered, non-steady state intermediates
		return production_rates[self.scheme._mask_stable]

	#Simulate the kinetic system up to a t_max
	def simulate(self,t_max,N_grid=100):
		t = np.linspace(0, t_max, N_grid)
		sol = odeint(self._sode, self.scheme._c_initial_stable, t)
		
		c_buffered = {
			i: {
				'x': t, 
				'y': np.array([self.scheme._c_initial_buffered[cnt]] * N_grid)
			} 
			for cnt, i in enumerate(self.scheme._labels_buffered)
		}

		c_intermediates = np.array([self._get_intermediates(i) for i in sol])
		c_intermediates = {i: {'x': t, 'y': j} for i,j in zip(self.scheme._labels_intermediates,c_intermediates.transpose())}
		c_stable = {i: {'x': t, 'y': j} for i,j in zip(self.scheme._labels_stable,sol.transpose())}
		
		#Save the result
		self.result = c_buffered | c_intermediates | c_stable
		
		print("Steady-state batch reactor simulation finished sucessfully")

class PLUG_FLOW_REACTOR_QSS():
	def __init__(self, scheme):
		self.scheme = scheme

	# same as for batch reactor
	def _soe(self,c_intermediates,c_stable):
		c = np.zeros(len(self.scheme._species))
		c[self.scheme._mask_buffered] = self.scheme._c_initial_buffered
		c[self.scheme._mask_intermediates] = c_intermediates
		c[self.scheme._mask_stable] = c_stable

		production_rates = self.scheme._get_all_prod_rate(c)
		return production_rates[self.scheme._mask_intermediates]

	# same as for batch reactor
	def _get_intermediates(self,c_stable):
		lgc_int_0 = np.array([-6]*self.scheme._N_intermediates)
		min_function = lambda lgc_int: np.sum([i**2 for i in self._soe(10**lgc_int,c_stable)])
		sol = minimize(min_function, lgc_int_0, method='Nelder-Mead', tol=1e-8, options={'maxiter': 1000, 'xatol': 1e-8, 'fatol': 1e-8, 'adaptive': True})

		if not sol.success:
			print("SOE for the reactive intermediates did not solve.")
			print(f"concentrations of stable species: {c_stable}")
			print(sol)
			return [np.Nan]*self.scheme._N_intermediates

		return 10**sol.x

	# differs from batch only by use of tau (( residence time = V/ v), v = volume flow rate))
	def _pfr_sode(self,c_stable, tau):
		c_intermediates = self._get_intermediates(c_stable)

		c = np.zeros(len(self.scheme._species))
		c[self.scheme._mask_buffered] = self.scheme._c_initial_buffered
		c[self.scheme._mask_intermediates] = c_intermediates
		c[self.scheme._mask_stable] = c_stable

		production_rates = self.scheme._get_all_prod_rate(c)
		return production_rates[self.scheme._mask_stable]

	# differs from batch only by use of tau 
	def simulate(self,tau_max,N_grid=100):
		tau = np.linspace(0, tau_max, N_grid)
		sol = odeint(self._pfr_sode, self.scheme._c_initial_stable, tau)
		
		c_buffered = {
			i: {
				'x': tau, 
				'y': np.array([self.scheme._c_initial_buffered[cnt]] * N_grid)
			} 
			for cnt, i in enumerate(self.scheme._labels_buffered)
		}

		c_intermediates = np.array([self._get_intermediates(i) for i in sol])
		c_intermediates = {i: {'x': tau, 'y': j} for i,j in zip(self.scheme._labels_intermediates,c_intermediates.transpose())}
		c_stable = {i: {'x': tau, 'y': j} for i,j in zip(self.scheme._labels_stable,sol.transpose())}
		
		#Save the result
		self.result = c_buffered | c_intermediates | c_stable

		print("Steady-state PFR simulation finished sucessfully")

class CSTR_REACTOR_QSS():
	def __init__(self, scheme):
		self.scheme = scheme

	# same as for batch reactor:
	def _soe(self,c_intermediates,c_stable):
		c = np.zeros(len(self.scheme._species))
		c[self.scheme._mask_buffered] = self.scheme._c_initial_buffered
		c[self.scheme._mask_intermediates] = c_intermediates
		c[self.scheme._mask_stable] = c_stable

		production_rates = self.scheme._get_all_prod_rate(c)
		return production_rates[self.scheme._mask_intermediates]
	
	def _get_intermediates(self,c_stable):
		lgc_int_0 = np.array([-6]*self.scheme._N_intermediates)
		min_function = lambda lgc_int: np.sum([i**2 for i in self._soe(10**lgc_int,c_stable)])
		sol = minimize(min_function, lgc_int_0, method='Nelder-Mead', tol=1e-8, options={'maxiter': 1000, 'xatol': 1e-8, 'fatol': 1e-8, 'adaptive': True})

		if not sol.success:
			print("SOE for the reactive intermediates did not solve.")
			print(f"concentrations of stable species: {c_stable}")
			print(sol)
			return [np.nan]*self.scheme._N_intermediates

		return 10**sol.x

	def _sode(self, c_stable, t, tau):
		c_intermediates = self._get_intermediates(c_stable)

		c = np.zeros(len(self.scheme._species))
		c[self.scheme._mask_buffered] = self.scheme._c_initial_buffered
		c[self.scheme._mask_intermediates] = c_intermediates
		c[self.scheme._mask_stable] = c_stable

		production_rates = self.scheme._get_all_prod_rate(c)

		# differs here - cstr rate = production rate + (C_0 - C)/tau
		flow_term = np.zeros(len(self.scheme._species))
		flow_term[self.scheme._mask_stable] =  (self.scheme._c_initial_stable - c_stable)/tau

		cstr_rates = production_rates[self.scheme._mask_stable] + flow_term[self.scheme._mask_stable]
		return cstr_rates

	def simulate(self, t_max, tau, N_grid=100):
		t = np.linspace(0, t_max, N_grid)

		sol = odeint(self._sode, self.scheme._c_initial_stable, t, args=(tau,))

		c_buffered = {
			i: {
				'x': t, 
				'y': np.array([self.scheme._c_initial_buffered[cnt]] * N_grid)
			} 
			for cnt, i in enumerate(self.scheme._labels_buffered)
		}

		c_intermediates = np.array([self._get_intermediates(i) for i in sol])
		c_intermediates = {i: {'x': t, 'y': j} for i,j in zip(self.scheme._labels_intermediates,c_intermediates.transpose())}
		c_stable = {i: {'x': t, 'y': j} for i,j in zip(self.scheme._labels_stable,sol.transpose())}
		
		self.result = c_buffered | c_intermediates | c_stable
		print("Steady-state CSTR simulation finished sucessfully")