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
a_CO2 = 2e3

#===================================================================================================================
#---------------------------------------------------------------------------------------Private methods
#===================================================================================================================
#-----------------------------------------------------------System of equations
def _soe(p,c0,c,T,a_CO2):
    x, y = p
    
    lgH2O = -40+42*x
    lgO2 = -120+122*y
        
    c['H2O'] = 10**lgH2O
    c['O2'] = 10**lgO2
    
    c['NO'] = c0['N']/(1+_reactions['NO/NO2']['K_p'](T)*c['O2']**0.5+_reactions['NO2/HNO2']['K_p'](T)*_reactions['NO/NO2']['K_p'](T)*c['O2']**0.5*c['H2O']**0.5/c['O2']**0.25+_reactions['HNO2/HNO3']['K_p'](T)*_reactions['NO2/HNO2']['K_p'](T)*_reactions['NO/NO2']['K_p'](T)*c['O2']**0.5*c['H2O']**0.5/c['O2']**0.25*c['O2']**0.5)
    c['NO2'] = _reactions['NO/NO2']['K_p'](T)*c['NO']*c['O2']**0.5
    c['HNO2'] = _reactions['NO2/HNO2']['K_p'](T)*c['NO2']*c['H2O']**0.5/c['O2']**0.25
    c['HNO3'] = _reactions['HNO2/HNO3']['K_p'](T)*c['HNO2']*c['O2']**0.5
    
    if c0['O']>=(2*c0['S']+c0['N']+c0['H']/2):
        c['SO2'] = c0['S']/(1+_reactions['SO2/SO3']['K_p'](T)*c['O2']**0.5+_reactions['SO3/H2SO4']['K_p'](T)*_reactions['SO2/SO3']['K_p'](T)*c['O2']**0.5*c['H2O'])
            
        c['SO3'] = Reactions['SO2/SO3']['K_p'](T)*c['SO2']*c['O2']**0.5
        c['H2SO4'] = Reactions['SO3/H2SO4']['K_p'](T)*c['SO3']*c['H2O']

        c['S'] = 0
        c['COS'] = 0
        c['H2S'] = 0
    else:
        c['SO2'] = _reactions['S/SO2']['K_p'](T)*c['O2']                      #<------------ This equation does not hold if there is no S!
        c['H2S'] = c['H2O']/c['O2']**0.5/_reactions['H2S/S']['K_p'](T)
        c['COS'] = _reactions['H2S/COS']['K_p'](T)*a_CO2*c['H2S']/c['H2O']    #<------------ This equation does not scale
        
        c['S'] = c0['S'] - (c['COS'] + c['H2S'] + c['SO2'])
        c['SO3'] = 0
        c['H2SO4'] = 0
        
    q = [2*c['H2O'] + 2*c['H2S'] + c['HNO2'] + c['HNO3'] + 2*c['H2SO4'] - c0['H'],
        2*c['O2'] + c['H2O'] - c['COS'] + c['NO'] + 2*c['NO2'] + 2*c['HNO2'] + 3*c['HNO3'] + 2*c['SO2'] + 3*c['SO3'] + 4*c['H2SO4'] - c0['O']]

    return np.array(q)
    
def _solve(c0,T=298.15,a_CO2=a_CO2,verbose=True,local=True):
    
    c = {prod: 0 for prod in _products}
    
    if c0['O']<c0['H']/2+c0['N']-c0['S']:                     #Outside of the studied range
        if verbose:
            print('Oxygen too low. Outside of studied range.')
        return c
    
    if local:
        x0 = np.linspace(0,1,10) 
        y0 = np.linspace(0,1,10)  

        points = [[x,y] for x in x0 for y in y0]               #1000 points or less do not increase the calculation time significantly
        energies = [np.sum(np.abs(_soe(p,c0,c,T,a_CO2))) for p in points]
        indx = np.argmin(energies)
        x0, y0 = points[indx]

        sol = minimize(lambda p: np.sum([i**2 for i in _soe(p,c0,c,T,a_CO2)]), [x0,y0], method='Nelder-Mead', tol=1e-6)
        #Reducing the tolerance will quickly lead to overflow errors
    else:
        sol = differential_evolution(lambda p: np.sum([i**2 for i in _soe(p,c0,c,T,a_CO2)]), [(0,1),(0,1)], tol=1e-12)
        #Increasing the tolerance will quickly lead to local solutions

    #Check if it solved
    if not sol.success:
        print(sol.message)
        return {prod: 0 for prod in _products}
        
    #Check if local minimum
    if np.all(np.abs(sol.fun)>1e-4):
        if local:
            if verbose:
                print('Convergence doubtful. Defaulting to global solve.')
            c = _solve(c0,T,a_CO2=a_CO2,verbose=verbose,local=False)
        else:
            if verbose:
                print('Convergence doubtful.')
    
    #-----------------------------------------------------------Output
    return c

def _solve_ppmx(p0,T=298.15,a_CO2=a_CO2,**kwargs):
    '''p0 = {'H2O': 10, 'H2S': 3, 'O2': 2, 'NO2': 2.5, 'SO2': 1} in ppmx'''
    
    c0 = {key: value*c_CO2*1e-6 for key,value in p0.items()}                                #[mol/m^3]

    c0['H'] = 2*c0['H2S']+2*c0['H2O']                                                       #[mol/m^3]
    c0['N'] = c0['NO2']                                                                     #[mol/m^3]
    c0['O'] = c0['H2O']+2*c0['SO2']+2*c0['NO2']+2*c0['O2']                                  #[mol/m^3]
    c0['S'] = c0['SO2']+c0['H2S']                                                           #[mol/m^3]
    
    c0 = {key: c0[key] for key in _elements}
    
    return _solve(c0,T,a_CO2,**kwargs)
    
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
        'T': temperature in [K]
    }'''
    
    if not isinstance(x0, dict):
        print('Wrong input!')
        return -1
        
    #Add default values, if not specified
    for key, default in [('CO2', 2e3),('T', 298.15)]:
        if not key in x0:
            x0 = x0 | {key: default}
        
    if not set(_elements) - set(x0.keys()):               #True if x has {'H','N','O','S'} keys
        for key in _elements:
            if not _domain[key]['min'] <= x0[key] <= _domain[key]['max']:
                print(f'Wrong input! {key} outside of range.')
                return -1
                
        sol = _solve(x0,T=x0['T'],a_CO2=x0['CO2'],**kwargs)
        
    elif not set(_impurities) - set(x0.keys()):           #True if x has {'H2O','H2S', 'O2','NO2','SO2'} keys
    
        for key in _impurities:
            if not _domain[key]['min'] <= x0[key] <= _domain[key]['max']:
                print(f'Wrong input! {key} outside of range.')
                return -1
                
        sol = _solve_ppmx(x0,T=x0['T'],a_CO2=x0['CO2'],**kwargs)
    else:
        print(f'Wrong input! Insufficient number of concentrations provided.')
        return -1
        
    return {key: float(c) for key,c in sol.items()}
    


