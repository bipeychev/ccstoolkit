#!/usr/bin/python3

#===================================================================================================================
#---------------------------------------------------------------------------------------Import libraries
#===================================================================================================================
import numpy as np
from scipy.optimize import minimize, differential_evolution

from . import _reactions as _reactions_

_domain = _reactions_.get_domain()
_reactions = _reactions_.get_reactions()

_elements = sorted(['H','N','O','S'])
_impurities = sorted(['H2O','H2S','O2','NO2','SO2'])
_products = sorted(['O2','H2O','COS','NO','NO2','HNO2','HNO3', 'H2S','S','SO2','SO3','H2SO4'])
_c_CO2 = 18.55e3
_a_CO2 = 2e3

#===================================================================================================================
#---------------------------------------------------------------------------------------Private methods
#===================================================================================================================
#-----------------------------------------------------------System of equations
def _soe(p,c0,c,T,a_CO2):
	'''System of ordinary equations (4 mass balances + 8 equilibrium conditions).'''
	#------------------Normalization/scaling of the minimizing parameters
	x, y = p
	
	lgH2O = -40+42*x
	lgO2 = -120+122*y
		
	c['H2O'] = 10**lgH2O
	c['O2'] = 10**lgO2
	
	#------------------Finding the concentration of nitrogen compounds
	c['NO'] = c0['N']/(1+_reactions['NO/NO2']['K_p'](T)*c['O2']**0.5+_reactions['NO2/HNO2']['K_p'](T)*_reactions['NO/NO2']['K_p'](T)*c['O2']**0.5*c['H2O']**0.5/c['O2']**0.25+_reactions['HNO2/HNO3']['K_p'](T)*_reactions['NO2/HNO2']['K_p'](T)*_reactions['NO/NO2']['K_p'](T)*c['O2']**0.5*c['H2O']**0.5/c['O2']**0.25*c['O2']**0.5)
	c['NO2'] = _reactions['NO/NO2']['K_p'](T)*c['NO']*c['O2']**0.5
	c['HNO2'] = _reactions['NO2/HNO2']['K_p'](T)*c['NO2']*c['H2O']**0.5/c['O2']**0.25
	c['HNO3'] = _reactions['HNO2/HNO3']['K_p'](T)*c['HNO2']*c['O2']**0.5
	
	#------------------Finding the concentration of sulphur compounds
	#---------Oxidative environment
	if c0['O']>=(2*c0['S']+c0['N']+c0['H']/2):
	
		c['SO2'] = c0['S']/(1+_reactions['SO2/SO3']['K_p'](T)*c['O2']**0.5+_reactions['SO3/H2SO4']['K_p'](T)*_reactions['SO2/SO3']['K_p'](T)*c['O2']**0.5*c['H2O'])
			
		c['SO3'] = _reactions['SO2/SO3']['K_p'](T)*c['SO2']*c['O2']**0.5
		c['H2SO4'] = _reactions['SO3/H2SO4']['K_p'](T)*c['SO3']*c['H2O']

		c['S'] = 0
		c['COS'] = 0
		c['H2S'] = 0
		
	#---------Reductive environment
	else:
	
		#This equation does not hold if there is no S!
		c['SO2'] = _reactions['S/SO2']['K_p'](T)*c['O2']
		#This equation does not scale
		c['H2S'] = c['H2O']/c['O2']**0.5/_reactions['H2S/S']['K_p'](T)
		c['COS'] = _reactions['H2S/COS']['K_p'](T)*a_CO2*c['H2S']/c['H2O']
		
		c['S'] = c0['S'] - (c['COS'] + c['H2S'] + c['SO2'])
		c['SO3'] = 0
		c['H2SO4'] = 0
	
	#------------------Remaining mass balances {H,O}
	q = [2*c['H2O'] + 2*c['H2S'] + c['HNO2'] + c['HNO3'] + 2*c['H2SO4'] - c0['H'],
		2*c['O2'] + c['H2O'] - c['COS'] + c['NO'] + 2*c['NO2'] + 2*c['HNO2'] + 3*c['HNO3'] + 2*c['SO2'] + 3*c['SO3'] + 4*c['H2SO4'] - c0['O']]

	#------------------Return
	return np.array(q)
	
def _solve(c0,T=298.15,a_CO2=_a_CO2,verbose=True,local=True,**kwargs):
	'''Find the equilibrium composition given the total concentration of {H,N,O,S}, temperature, and CO2 activity.'''
	
	#------------------Initialize the dictionary for the results
	c = {prod: 0 for prod in _products}
	
	#------------------Check if inside the studied range
	if c0['O']<c0['H']/2+c0['N']-c0['S']:					 #Outside of the studied range
		if verbose:
			print('Oxygen too low. Outside of studied range.')
		return -1
	
	#------------------Local solve
	if local:
	
		#------------------Find rough initial guess
		x0 = np.linspace(0,1,10) 
		y0 = np.linspace(0,1,10)  

		#1000 points or less do not increase the calculation time significantly
		points = [[x,y] for x in x0 for y in y0]			   
		energies = [np.sum(np.abs(_soe(p,c0,c,T,a_CO2))) for p in points]
		indx = np.argmin(energies)
		x0, y0 = points[indx]
		
		#------------------Find local solution
		#Reducing the tolerance will quickly lead to overflow errors
		sol = minimize(lambda p: np.sum([i**2 for i in _soe(p,c0,c,T,a_CO2)]), [x0,y0], method='Nelder-Mead', tol=1e-6, **kwargs)
		
	#------------------Global solve
	else:
	
		#Increasing the tolerance will quickly lead to local solutions
		sol = differential_evolution(lambda p: np.sum([i**2 for i in _soe(p,c0,c,T,a_CO2)]), [(0,1),(0,1)], tol=1e-12, **kwargs)

	#------------------Check if it solved
	if not sol.success:
		print(sol.message)
		return {prod: 0 for prod in _products}
		
	#------------------Check if local minimum
	if np.all(np.abs(sol.fun)>1e-4):			 #Could give a false positive?
		if local:
			if verbose:
				print('Convergence doubtful. Defaulting to global solve.')
			c = _solve(c0,T=T,a_CO2=a_CO2,verbose=verbose,local=False)
		else:
			if verbose:
				print('Convergence doubtful.')
	
	#------------------Output
	return c

def _solve_ppmx(p0,T=298.15,a_CO2=_a_CO2,c_tot=_c_CO2,**kwargs):
	'''Find the equilibrium composition given the initial concentration of {H2O,H2S,O2,NO2,SO2} in ppm, temperature, total concentration in mol/m^3, and CO2 activity.
	p0 = {'H2O': 10, 'H2S': 3, 'O2': 2, 'NO2': 2.5, 'SO2': 1} in ppmx'''
	
	#------------------Calculate elemental concentrations in mol/m^3
	c0 = _calc_element_conc_from_ppmx(p0,c_tot)
	
	#------------------Return the solution
	return _solve(c0,T=T,a_CO2=a_CO2,**kwargs)
	
def _convert_ppmx_to_mM(p0,c_tot=_c_CO2):
	return {key: value*c_tot*1e-6 for key,value in p0.items()}							  #[mol/m^3]
	
def _calc_element_conc(c0):
	'''Calculate the concentration of {H,N,O,S} from the concentrations of {H2O,H2S,O2,NO2,SO2}.'''
	
	#------------------Initialize a dict for the results
	c = {key: 0 for key in _elements}
	
	#------------------Calculate the concentration of the elements
	c['H'] = 2*c0['H2S']+2*c0['H2O']													  #[mol/m^3]
	c['N'] = c0['NO2']																	  #[mol/m^3]
	c['O'] = c0['H2O']+2*c0['SO2']+2*c0['NO2']+2*c0['O2']								  #[mol/m^3]
	c['S'] = c0['SO2']+c0['H2S']														  #[mol/m^3]
	
	#------------------Return
	return c
	
def _calc_element_conc_from_ppmx(p0,c_tot=_c_CO2):

	#------------------Convert to mol/m^3
	c0 = _convert_ppmx_to_mM(p0,c_tot)
	
	#------------------Calculate elemental concentrations
	c0 = _calc_element_conc(c0)
	
	#------------------Return
	return c0
	
#===================================================================================================================
#---------------------------------------------------------------------------------------Public methods
#===================================================================================================================
def get_composition(x0,**kwargs):
	'''x0 = {
		'H': total hydrogen concentration in [mM], 
		'N': total nitrogen concentration in [mM], 
		'O': total oxygen concentration in [mM], 
		'S': total sulphur concentration in [mM],
		'CO2': activity of CO2 in [mM],
		'T': temperature in [K]
	}
	
	x0 = {
		'H2O': concentration in [ppmx], 
		'H2S': concentration in [ppmx], 
		'O2': concentration in [ppmx], 
		'NO2': concentration in [ppmx], 
		'SO2': concentration in [ppmx],
		'CO2': activity of CO2 in [mM],
		'tot': total concentration of all species [mM],
		'T': temperature in [K]
	}'''
	
	#------------------Check if the input is a dict
	if not isinstance(x0, dict):
		print('Wrong input!')
		return -1
		
	#------------------Add default values, if not specified
	for key, default in [('CO2', _a_CO2),('tot', _c_CO2),('T', 298.15)]:
		if key not in x0:
			x0 = x0 | {key: default}
	
	#------------------Solve given {H,N,O,S}	
	#True if x0 has {'H','N','O','S'} keys
	if not set(_elements) - set(x0.keys()):	
	
		#------------------Check if the concentrations are within the specified range
		for key in _elements:
			if not _domain[key]['min'] <= x0[key] <= _domain[key]['max']:
				print(f'Wrong input! {key} outside of range.')
				return -1
				
		#------------------Solve		
		sol = _solve(x0,T=x0['T'],a_CO2=x0['CO2'],**kwargs)
	
	#------------------Solve given {H2O,H2S,O2,NO2,SO2}	
	#True if x0 has any of {'H2O','H2S', 'O2','NO2','SO2'}
	elif set(_impurities) & x0.keys():  #not set(_impurities) - set(x0.keys()):
	
		#------------------Add default values, if not specified
		for key in _impurities:
			if key not in x0:
				x0 = x0 | {key: 0}
	
		#------------------Calculate element concentrations
		c0 = _calc_element_conc_from_ppmx(x0,x0['tot'])
		
		#------------------Check if the concentrations are within the specified range
		for key in _elements:
			if not _domain[key]['min'] <= c0[key] <= _domain[key]['max']:
				print(f'Wrong input! {key} outside of range.')
				return -1
				
		#------------------Solve			
		sol = _solve(c0,T=x0['T'],a_CO2=x0['CO2'],**kwargs)
		
	#------------------Wrong input
	else:
		print('Wrong input! Insufficient number of concentrations provided.')
		return -1
	
	#------------------Return the equilibrium composition	
	#True if the oxygen is too low
	if sol == -1:
		return -2
	else:
		return {key: float(c) for key,c in sol.items()}
	
