#!/usr/bin/python3

#===================================================================================================================
#---------------------------------------------------------------------------------------Import
#===================================================================================================================
import re
import numpy as np

import ccstoolkit.common._line_logic as _line_logic
import ccstoolkit.common._substances as _substances_

_substances = _substances_.get_substances_TD_data()

R = 8.314                                                                                   #[J/K/mol]

#===================================================================================================================
#---------------------------------------------------------------------------------------Data
#===================================================================================================================
#'...at room T, the saturation point of water is 1.3 mM (3.2 kPa)...'
#'On the other hand, meta stable supersaturated mixtures are possible, 
#so we plot the diagrams to an upper limit of lgH2O = 2.'
_bounds = {
    'x': (-20, 2),                                                                          #[-]
    'y': (-120, 2)                                                                          #[-]
}

_domain = {
    'H': {'min': 0.015, 'max': 12},                                                     #[mol/m^3]
    'N': {'min': 0.015, 'max': 4},                                                      #[mol/m^3]
    'O': {'min': 0.015, 'max':24},                                                      #[mol/m^3]
    'S': {'min': 0.015, 'max':4},                                                       #[mol/m^3]
    'CO2': {'min': 0.015, 'max':3e3},                                                   #[mol/m^3]
    'H2O': {'min': 1, 'max':200},                                                       #[ppmx]
    'H2S': {'min': 1, 'max':100},                                                       #[ppmx]
    'NO2': {'min': 1, 'max':200},                                                       #[ppmx]
    'O2': {'min': 1, 'max':200},                                                        #[ppmx]
    'SO2': {'min': 1, 'max':100},                                                       #[ppmx]
    'T': {'min': 273.15-50, 'max': 273.15+100}                                          #[K]
}

_reactions = {
    'H2S/S': {'reaction': {'reaction': 'H2S + 1/2O2  <=> S + H2O'}},
    'H2S/COS': {'reaction': {'reaction': 'CO2 + H2S  <=> COS + H2O'}},
    'S/SO2': {'reaction': {'reaction': 'S + O2  <=> SO2'}},
    'SO2/SO3': {'reaction': {'reaction': 'SO2 + 1/2O2  <=> SO3'}},
    'SO3/H2SO4': {'reaction': {'reaction': 'SO3 + H2O  <=> H2SO4'}},
    'NO/NO2': {'reaction': {'reaction': 'NO + 1/2O2 <=> NO2'}},
    'NO2/HNO2': {'reaction': {'reaction': 'NO2 + 1/2H2O <=> HNO2 + 1/4O2'}},
    'HNO2/HNO3': {'reaction': {'reaction': 'HNO2 + 1/2O2 <=> HNO3'}},
    'COS/S': {'reaction': {'reaction': 'CO2 + S <=> COS + 1/2O2'}},
    'SO2/H2SO4': {'reaction': {'reaction': 'SO2 + H2O + 1/2O2 <=> H2SO4'}},
    'NO/HNO2': {'reaction': {'reaction': 'NO + 1/2H2O + 1/4O2 <=> HNO2'}},
    'NO2/HNO3': {'reaction': {'reaction': 'NO2 + 1/2H2O + 1/4O2 <=> HNO3'}},
    'H2S/H2SO4': {'reaction': {'reaction': 'H2S + 2O2 <=> H2SO4'}},
    'S/H2SO4': {'reaction': {'reaction': 'S + H2O + 3/2O2 <=> H2SO4'}},
    'SO2/H2S': {'reaction': {'reaction': 'SO2 + H2O <=> H2S + 3/2O2'}}
}

_rules = {
    'H2S/S+NO': [('right of','H2S/COS+NO'),('below','S/SO2+NO'),('below','H2S/H2SO4+NO')],                         
    'H2S/COS+NO': [('below','COS/S+NO')],
    'S/SO2+NO': [('left of','SO2/H2SO4+NO'),('left of','H2S/S+NO')],
    'SO2/SO3+NO': [('left of','SO3/H2SO4+NO')],
    'H2S/H2SO4+NO': [('right of','SO2/H2SO4+NO'),('right of','H2S/S+NO')],
    'COS/S+NO': [('left of','H2S/COS+NO')],
    'SO2/H2SO4+NO': [('below', 'SO2/SO3+NO'),('above','S/SO2+NO'),('above','H2S/H2SO4+NO')],                
    'S/H2SO4+NO': [('below','S/SO2+NO'),('above','H2S/H2SO4+NO')],
    'SO2/H2S+NO': [('above','S/SO2+NO'),('below','H2S/H2SO4+NO')],
    'SO3/H2SO4+NO': [('above','SO2/SO3+NO'),('below','NO/NO2+SO3')],
    'SO3/H2SO4+NO2': [('above','SO2/SO3+NO'),('above','NO/NO2+SO3')],
    
    'NO/HNO2+H2SO4': [('below','NO/NO2+SO3')],
    'NO2/HNO3+H2SO4': [('above','HNO2/HNO3+H2SO4')],
    'NO2/HNO2+H2SO4': [('above','NO/HNO2+H2SO4'),('below','HNO2/HNO3+H2SO4')],
    'HNO2/HNO3+H2SO4': [('below','NO2/HNO2+H2SO4')],
    'NO/NO2+SO3': [('below','NO/HNO2+H2SO4'),('left of','SO3/H2SO4+NO')],
    'NO/NO2+H2SO4': [('below','NO/HNO2+H2SO4'),('right of','SO3/H2SO4+NO')]
}

#===================================================================================================================
#---------------------------------------------------------------------------------------Calculate/format/etc.
#===================================================================================================================
#Parse the reactions
for key, reaction in _reactions.items():
    reaction['key'] = key
    
    #---------------------------Extract the reactants, product and coefficients---------------------------
    r = reaction['reaction']['reaction']
    r = r.replace(' ','').split('<=>')
    lhs, rhs = r
    
    #Reactants
    #'^(\d+(?:/\d+)?)?([A-Za-z].*)' matches all digits before the first non-digit (except '/')
    lhs_coeffs = [re.match(r'^(\d+(?:/\d+)?)?([A-Za-z].*)', i).group(1) for i in lhs.split('+')]
    lhs_coeffs = [(int(coeff.split('/')[0])/int(coeff.split('/')[1]) if '/' in coeff else int(coeff)) if coeff else 1 for coeff in lhs_coeffs]
    lhs_substances = [re.match(r'^(\d+(?:/\d+)?)?([A-Za-z].*)', i).group(2) for i in lhs.split('+')]
    
    reaction['reaction']['lhs'] = {'substances': lhs_substances, 'coeffs': lhs_coeffs}
    
    #Products
    rhs_coeffs = [re.match(r'^(\d+(?:/\d+)?)?([A-Za-z].*)', i).group(1) for i in rhs.split('+')]
    rhs_coeffs = [(int(coeff.split('/')[0])/int(coeff.split('/')[1]) if '/' in coeff else int(coeff)) if coeff else 1 for coeff in rhs_coeffs]
    rhs_substances = [re.match(r'^(\d+(?:/\d+)?)?([A-Za-z].*)', i).group(2) for i in rhs.split('+')]
    
    reaction['reaction']['rhs'] = {'substances': rhs_substances, 'coeffs': rhs_coeffs}
    
    #All coefficients (with a minus for the reactants) and substances
    Coeffs = [-i for i in lhs_coeffs]+rhs_coeffs
    Subs = lhs_substances+rhs_substances
    reaction['reaction']['coeffs'] = Coeffs
    reaction['reaction']['substances'] = Subs
    
    #---------------------------Calculate the free energy of the reaction---------------------------
    drg = 0
    drh = 0
    drcp = 0
    nur = 0
    for s,c in zip(Subs,Coeffs):
        drg += _substances[s]['dfg']*c
        drh += _substances[s]['dfh']*c
        drcp += _substances[s]['cp']*c
        nur += c*(not _substances[s]['solid'])
        
    reaction['drg'] = drg*1e3                               #[J/mol]
    reaction['drh'] = drh*1e3                               #[J/mol]
    reaction['drcp'] = drcp                                 #[J/mol/K]
    reaction['nur'] = nur
    
    #---------------------------Calculate the equilibrium constant of the reaction---------------------------
    reaction['K_chi_298'] = np.exp(-(reaction['drg'])/R/298.15)
    reaction['K_chi'] = lambda T, reaction=reaction: reaction['K_chi_298']*(T/298.15)**(reaction['drcp']/R)*np.exp(-(reaction['drh']-298.15*reaction['drcp'])/R*(298.15-T)/T/298.15)
    
    reaction['K_p'] = lambda T, reaction=reaction: reaction['K_chi'](T)*(1e5/R/T)**reaction['nur']
    reaction['lgK_p'] = lambda T, reaction=reaction: np.log10(reaction['K_p'](T))
    
    #---------------------------Equation in the from a+b*lgH2O+c*lgO2=0---------------------------
    ind = Subs.index('O2') if 'O2' in Subs else None
    coeff_O2 = Coeffs[ind] if ind!=None else 0

    ind = Subs.index('H2O') if 'H2O' in Subs else None
    coeff_H2O = Coeffs[ind] if ind!=None else 0

    ind = Subs.index('CO2') if 'CO2' in Subs else None
    coeff_CO2 = Coeffs[ind] if ind!=None else 0
    
    #Substances containing N
    inds = [Subs.index(i) if i in Subs else None for i in ['NO','NO2','HNO2','HNO3']]
    #All of them are 1/2 of N_tot (there are 2 N containing substances per reaction that are equally stable)
    coeff_Ntot = sum([Coeffs[ind] if ind!=None else 0 for ind in inds])
    
    #Substances containing N
    inds = [Subs.index(i) if i in Subs else None for i in ['COS','H2S','SO2','SO3','H2SO4']]   #S is solid!
    #All of them are 1/2 of S_tot (there are 2 S containing substances per reaction that are equally stable)
    coeff_Stot = sum([Coeffs[ind] if ind!=None else 0 for ind in inds])
    
    #P = {'S': 1, 'N': 1, 'C': 2000, 'T': 298.15}   <---------   Important!!!
    #The intercept of the equations is a function of the composition.
    intercept = lambda P, lgKp=_reactions[key]['lgK_p'], c_s=coeff_Stot, c_n=coeff_Ntot, c_c=coeff_CO2: -lgKp(P['T'])+c_s*np.log10(1/2*P['S'])+c_n*np.log10(1/2*P['N'])+c_c*np.log10(P['CO2'])
    slope_x = coeff_H2O
    slope_y = coeff_O2
    
    reaction['line'] = {'equation': 'c[0](P)+c[1]*lgH2O+c[2]*lgO2=0 '}
    reaction['line']['coeffs'] = lambda P, c=(intercept,slope_x,slope_y): (c[0](P),c[1],c[2])
    
    #---------------------------Equations in the form lgH2O=f(lgO2) and lgO2=f(lgH2O)---------------------------
    reaction['line']['x'] = lambda y, P, c=(intercept,slope_x,slope_y): -(c[0](P)+c[2]*y)/c[1]
    reaction['line']['y'] = lambda x, P, c=(intercept,slope_x,slope_y): -(c[0](P)+c[1]*x)/c[2]
    
    #---------------------------Is the line horizontal/vertical---------------------------
    reaction['line']['vertical'] = coeff_O2 == 0
    reaction['line']['horizontal'] = coeff_H2O == 0 
   
#Create the lines dicts
_unformatted_lines = {key: _reactions[key.split('+')[0]]['line'] for key,item in _rules.items()}
_lines = _line_logic._form_lines(_rules,_unformatted_lines)

#===================================================================================================================
#---------------------------------------------------------------------------------------Export
#===================================================================================================================
def get_domain():
	return _domain
	
def get_reactions():
	return _reactions

def get_bounds():
    return _bounds
    
def get_lines():
    return _lines


