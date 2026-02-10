#!/usr/bin/python3

#===================================================================================================================
#---------------------------------------------------------------------------------------Import
#===================================================================================================================
import re
import numpy as np

from . import _math
from . import _substances

_substances = _substances.get_substances_TD_data()

R = 8.314                                                                                   #[J/K/mol]

#===================================================================================================================
#---------------------------------------------------------------------------------------Data
#===================================================================================================================
#'...at room T, the saturation point of water is 1.3 mM (3.2 kPa)...'
#'On the other hand, meta stable supersaturated mixtures are possible, 
#so we plot the diagrams to an upper limit of lgH2O = 2.'
_bounds = {
    'x': (-20, 2),
    'y': (-120, 2)
}

_domain = {
    'S': {'min': 0.015, 'max': 3},
    'N': {'min': 0.015, 'max': 3},
    'C': {'min': 0.015, 'max':3e3},
    'T': {'min': 258.15, 'max': 373.15}
}

#Specify the reactions
_reactions = {
    'Fe/Fe3O4': {'reaction': {'reaction': '3Fe + 2O2 <=> Fe3O4'}},
    'Fe/Fe(OH)2': {'reaction': {'reaction': 'Fe + H2O + 1/2O2 <=> Fe(OH)2'}},
    'Fe3O4/Fe2O3': {'reaction': {'reaction': '2Fe3O4 + 1/2O2 <=> 3Fe2O3'}},
    'Fe3O4/Fe(OH)2': {'reaction': {'reaction': 'Fe3O4 + 3H2O <=> 3Fe(OH)2 + 1/2O2'}},
    'Fe2O3/FeO(OH)': {'reaction': {'reaction': 'Fe2O3 + H2O <=> 2FeO(OH)'}},
    'Fe3O4/FeO(OH)': {'reaction': {'reaction': 'Fe3O4 + 3/2H2O + 1/4O2 <=> 3FeO(OH)'}},
    'FeCO3/FeO(OH)': {'reaction': {'reaction': 'FeCO3 + 1/2H2O + 1/4O2 <=> FeO(OH) + CO2'}},
    'Fe/CO2/FeCO3': {'reaction': {'reaction': 'Fe + CO2 + 1/2O2 <=> FeCO3'}},
    'Fe/CO/FeCO3': {'reaction': {'reaction': 'Fe + CO + O2 <=> FeCO3'}},
    'FeCO3/Fe2O3': {'reaction': {'reaction': '2FeCO3 + 1/2O2 <=> Fe2O3 + 2CO2'}},
    'Fe2O3/SO3/Fe2(SO4)3': {'reaction': {'reaction': 'Fe2O3 + 3SO3 <=> Fe2(SO4)3'}},
    'Fe/COS/FeS': {'reaction': {'reaction': 'Fe + COS + 1/2O2 <=> FeS + CO2'}},
    'Fe/H2S/FeS': {'reaction': {'reaction': 'Fe + H2S + 1/2O2 <=> FeS + H2O'}},
    'FeCO3/H2S/FeS': {'reaction': {'reaction': 'FeCO3 + H2S <=> FeS + H2O + CO2'}},
    'FeS/S/FeCO3': {'reaction': {'reaction': 'FeS + CO2 + 1/2O2 <=> FeCO3 + S'}},
    'FeCO3/S/FeSO4': {'reaction': {'reaction': 'FeCO3 + S + 3/2O2 <=> FeSO4 + CO2'}},
    'FeS/COS/FeS2': {'reaction': {'reaction': 'FeS + COS + 1/2O2 <=> FeS2 + CO2'}},
    'FeS/H2S/FeS2': {'reaction': {'reaction': 'FeS + H2S + 1/2O2 <=> FeS2 + H2O'}},
    'FeS2/S/FeSO4': {'reaction': {'reaction': 'FeS2 + 2O2 <=> FeSO4 + S'}},
    'FeSO4/S/Fe2(SO4)3': {'reaction': {'reaction': '2FeSO4 + S + 2O2 <=> Fe2(SO4)3'}},
    'FeSO4/SO2/Fe2(SO4)3': {'reaction': {'reaction': '2FeSO4 + SO2 + O2 <=> Fe2(SO4)3'}},
    'Fe2O3/NO2/Fe(NO3)2': {'reaction': {'reaction': 'Fe2O3 + 4NO2 + 1/2O2 <=> 2Fe(NO3)2'}},
    'Fe(NO3)2/NO2/FeO(OH)': {'reaction': {'reaction': 'Fe(NO3)2 + 1/2H2O <=> FeO(OH) + 2NO2 + 1/4O2'}},
    'Fe(NO3)2/HNO3/FeO(OH)': {'reaction': {'reaction': 'Fe(NO3)2 + 3/2H2O + 1/4O2 <=> FeO(OH) + 2HNO3'}},
    'FeSO4.H2O/FeSO4': {'reaction': {'reaction': 'FeSO4.H2O <=> FeSO4 + H2O'}},
    'FeSO4.7H2O/FeSO4.H2O': {'reaction': {'reaction': 'FeSO4.7H2O <=> FeSO4.H2O + 6H2O'}},
    'FeSO4.H2O/FeS2': {'reaction': {'reaction': 'FeSO4.H2O + S <=> FeS2 + 2O2 + H2O'}},
    'FeSO4.H2O/S/Fe2(SO4)3': {'reaction': {'reaction': '2FeSO4.H2O + S + 2O2 <=> Fe2(SO4)3 + 2H2O'}},
    'FeSO4.7H2O/S/FeS2': {'reaction': {'reaction': 'FeSO4.7H2O + S <=> FeS2 + 2O2 + 7H2O'}},
    'FeSO4.H2O/SO2/Fe2(SO4)3': {'reaction': {'reaction': '2FeSO4.H2O + SO2 + O2 <=> Fe2(SO4)3 + 2H2O'}},
    'FeSO4.7H2O/SO2/Fe2(SO4)3': {'reaction': {'reaction': '2FeSO4.7H2O + SO2 + O2 <=> Fe2(SO4)3 + 14H2O'}},
    'FeSO4.7H2O/H2SO4/Fe2(SO4)3': {'reaction': {'reaction': '2FeSO4.7H2O + H2SO4 + 1/2O2 <=> Fe2(SO4)3 + 15H2O'}},
    'FeSO4.7H2O/H2S/FeS2': {'reaction': {'reaction': 'FeSO4.7H2O + H2S  <=> FeS2 + 8H2O + 3/2O2'}},
    'FeSO4.7H2O/FeS': {'reaction': {'reaction': 'FeSO4.7H2O  <=> FeS + 7H2O + 2O2'}},
    'FeSO4.7H2O/Fe': {'reaction': {'reaction': 'FeSO4.7H2O  <=> Fe + H2S + 6H2O + 5/2O2'}},
    'FeCO3/Fe3O4': {'reaction': {'reaction': '3FeCO3 + 1/2O2 <=> Fe3O4 + 3CO2'}},
    'FeCO3/Fe(OH)2': {'reaction': {'reaction': 'FeCO3 + H2O <=> Fe(OH)2 + CO2'}},
    'Fe(OH)2/FeO(OH)': {'reaction': {'reaction': 'Fe(OH)2 + 1/4O2 <=> FeO(OH) + 1/2H2O'}},
    'Fe(NO3)2/HNO2/FeO(OH)': {'reaction': {'reaction': 'Fe(NO3)2 + 3/2H2O  <=> FeO(OH) + 2HNO2 + 3/4O2'}},
    'Fe(NO3)2/NO/FeO(OH)': {'reaction': {'reaction': 'Fe(NO3)2 + 1/2H2O <=> FeO(OH) + 2NO + 5/4O2'}},
    'Fe2O3/NO/Fe(NO3)2': {'reaction': {'reaction': 'Fe2O3 + 4NO + 5/2O2 <=> 2Fe(NO3)2'}},
    'Fe2O3/HNO3/Fe(NO3)2': {'reaction': {'reaction': 'Fe2O3 + 4HNO3 <=> 2Fe(NO3)2 + 2H2O + 1/2O2'}}
}

#Specify each graph as a set of rules bounding the domain of each line
#Only lines with a rule will be included in the graph (the rules can be an empty list)
_rules_Fe_O = {                                 #Corrosion products of steel/Fe in the presence of oxygen
    'Fe/Fe3O4': [('below','Fe/Fe(OH)2')],
    'Fe/Fe(OH)2': [('below','Fe/Fe3O4')],
    'Fe3O4/Fe(OH)2': [('above','Fe/Fe3O4'),('below','Fe(OH)2/FeO(OH)')],
    'Fe3O4/Fe2O3': [('left of','Fe2O3/FeO(OH)')],
    'Fe3O4/FeO(OH)': [('below','Fe3O4/Fe2O3'),('above','Fe3O4/Fe(OH)2')],
    'Fe2O3/FeO(OH)': [('above','Fe3O4/Fe2O3')],
    'Fe(OH)2/FeO(OH)': [('below','Fe3O4/Fe(OH)2')]
}

_rules_Fe_C = {                                  #Corrosion products of steel/Fe in the presence of oxygen and CO2
    'Fe/Fe3O4': [('below','Fe/CO2/FeCO3'),('below','Fe/Fe(OH)2')],
    'Fe/Fe(OH)2': [('right of', 'FeCO3/Fe(OH)2'),('below', 'Fe/Fe3O4')],
    'Fe3O4/Fe(OH)2': [('right of', 'FeCO3/Fe(OH)2'),('below','Fe(OH)2/FeO(OH)'),('above','Fe/Fe(OH)2')],
    'Fe3O4/Fe2O3': [('left of','Fe2O3/FeO(OH)'),('above','FeCO3/Fe3O4')],
    'Fe3O4/FeO(OH)': [('below','Fe3O4/Fe2O3'),('above','Fe3O4/Fe(OH)2'),('above','FeCO3/Fe3O4')],
    'Fe2O3/FeO(OH)': [('above','Fe3O4/Fe2O3'),('above','FeCO3/Fe2O3')],
    'Fe/CO2/FeCO3': [('left of', 'FeCO3/Fe(OH)2'),('below','Fe/Fe3O4')],
    'FeCO3/Fe2O3': [('above','Fe3O4/Fe2O3'),('left of', 'Fe2O3/FeO(OH)')],
    'FeCO3/FeO(OH)': [('right of', 'Fe2O3/FeO(OH)'),('left of', 'FeCO3/Fe(OH)2'),('below', 'FeCO3/Fe3O4')],
    'FeCO3/Fe3O4': [('left of', 'FeCO3/Fe(OH)2'),('below', 'Fe3O4/Fe2O3'),('below', 'Fe3O4/FeO(OH)'),('above', 'Fe/CO2/FeCO3')],
    'FeCO3/Fe(OH)2': [('below','FeCO3/Fe3O4'),('above','Fe/CO2/FeCO3'),('below','FeCO3/FeO(OH)')],
    'Fe(OH)2/FeO(OH)': [('below','Fe3O4/Fe(OH)2'),('right of', 'FeCO3/Fe(OH)2')]
}

_rules_Fe_N = {                                  #Corrosion products of steel/Fe in the presence of oxygen and N oxidants
    'Fe/Fe3O4': [('below','Fe/CO2/FeCO3'),('below','Fe/Fe(OH)2')],
    'Fe/Fe(OH)2': [('right of', 'FeCO3/Fe(OH)2'),('below', 'Fe/Fe3O4')],
    'Fe3O4/Fe(OH)2': [('right of', 'FeCO3/Fe(OH)2'),('below','Fe(OH)2/FeO(OH)'),('above','Fe/Fe(OH)2')],
    'Fe3O4/Fe2O3': [('left of','Fe2O3/FeO(OH)'),('above','FeCO3/Fe3O4')],
    'Fe3O4/FeO(OH)': [('below','Fe3O4/Fe2O3'),('above','Fe3O4/Fe(OH)2'),('above','FeCO3/Fe3O4')],
    'Fe2O3/FeO(OH)': [('above','Fe3O4/Fe2O3'),('above','FeCO3/Fe2O3'),('special')],
    'Fe/CO2/FeCO3': [('left of', 'FeCO3/Fe(OH)2'),('below','Fe/Fe3O4')],
    'FeCO3/Fe2O3': [('above','Fe3O4/Fe2O3'),('left of', 'Fe2O3/FeO(OH)')],
    'FeCO3/FeO(OH)': [('right of', 'Fe2O3/FeO(OH)'),('left of', 'FeCO3/Fe(OH)2'),('below', 'FeCO3/Fe3O4')],
    'FeCO3/Fe3O4': [('left of', 'FeCO3/Fe(OH)2'),('below', 'Fe3O4/Fe2O3'),('below', 'Fe3O4/FeO(OH)'),('above', 'Fe/CO2/FeCO3')],
    'FeCO3/Fe(OH)2': [('below','FeCO3/Fe3O4'),('above','Fe/CO2/FeCO3'),('below','FeCO3/FeO(OH)')],
    'Fe(OH)2/FeO(OH)': [('below','Fe3O4/Fe(OH)2'),('right of', 'FeCO3/Fe(OH)2')],
    'Fe2O3/NO2/Fe(NO3)2': [('left of','Fe2O3/FeO(OH)'),('above','Fe2O3/NO/Fe(NO3)2'),
                           ('below','Fe(NO3)2/HNO3/FeO(OH)'),('below','Fe2O3/HNO3/Fe(NO3)2')],
    'Fe(NO3)2/NO2/FeO(OH)': [('above','Fe2O3/NO2/Fe(NO3)2'),('below','Fe(NO3)2/HNO3/FeO(OH)'),('above','Fe(NO3)2/NO/FeO(OH)')],
    'Fe(NO3)2/HNO3/FeO(OH)': [('above','Fe(NO3)2/NO2/FeO(OH)'),('above','Fe2O3/NO2/Fe(NO3)2'),('special')],
    #'Fe(NO3)2/HNO2/FeO(OH)': [],
    'Fe(NO3)2/NO/FeO(OH)': [('above','Fe2O3/NO/Fe(NO3)2'),('above','Fe(NO3)2/NO2/FeO(OH)')],
    'Fe2O3/NO/Fe(NO3)2': [('above','Fe(NO3)2/NO/FeO(OH)'),('above','Fe2O3/NO2/Fe(NO3)2')],
    'Fe2O3/HNO3/Fe(NO3)2': [('above','Fe2O3/NO2/Fe(NO3)2'),('special')]
}

_rules_Fe_S = {                                  #Corrosion products of steel/Fe in the presence of oxygen and S oxidants
    'Fe/COS/FeS': [('above', 'Fe/H2S/FeS')],
    'Fe/H2S/FeS': [('above', 'Fe/COS/FeS')],
    'FeS/COS/FeS2': [('above', 'FeS/H2S/FeS2')],
    'FeS/H2S/FeS2': [('above', 'FeS/COS/FeS2')],
    'FeS2/S/FeSO4': [('left of', 'FeSO4.H2O/FeSO4')],
    'FeSO4/S/Fe2(SO4)3': [('left of', 'FeSO4.H2O/FeSO4'),('above', 'FeSO4/SO2/Fe2(SO4)3')],    
    'FeSO4/SO2/Fe2(SO4)3': [('left of', 'FeSO4.H2O/FeSO4'),('above', 'FeSO4/S/Fe2(SO4)3')],
    'FeSO4.H2O/FeSO4': [('above', 'FeS2/S/FeSO4'), ('special')],
    'FeSO4.7H2O/FeSO4.H2O': [('above', 'FeSO4.H2O/FeS2'), ('below', 'FeSO4.H2O/SO2/Fe2(SO4)3')],
    'FeSO4.H2O/FeS2': [('right of', 'FeSO4.H2O/FeSO4'), ('left of', 'FeSO4.7H2O/FeSO4.H2O')],
    'FeSO4.H2O/S/Fe2(SO4)3': [('right of', 'FeSO4.H2O/FeSO4'), ('above', 'FeSO4.H2O/SO2/Fe2(SO4)3'), ('left of', 'FeSO4.7H2O/FeSO4.H2O')],
    'FeSO4.7H2O/S/FeS2': [('right of', 'FeSO4.7H2O/FeSO4.H2O'), ('below', 'FeSO4.7H2O/H2S/FeS2')],
    'FeSO4.H2O/SO2/Fe2(SO4)3': [('right of', 'FeSO4.H2O/FeSO4'), ('left of', 'FeSO4.7H2O/FeSO4.H2O'), ('above', 'FeSO4.H2O/S/Fe2(SO4)3')],
    'FeSO4.7H2O/SO2/Fe2(SO4)3': [('above', 'FeSO4.H2O/SO2/Fe2(SO4)3'), ('above', 'FeSO4.H2O/S/Fe2(SO4)3'), ('above', 'FeSO4.7H2O/H2SO4/Fe2(SO4)3')],
    'FeSO4.7H2O/H2SO4/Fe2(SO4)3': [('above', 'FeSO4.7H2O/SO2/Fe2(SO4)3')],
    'FeSO4.7H2O/H2S/FeS2': [('below', 'FeSO4.7H2O/S/FeS2')]  
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
    coeffs = [-i for i in lhs_coeffs]+rhs_coeffs
    subs = lhs_substances+rhs_substances
    reaction['reaction']['coeffs'] = coeffs
    reaction['reaction']['substances'] = subs
    
    #---------------------------Calculate the free energy of the reaction---------------------------
    drg = 0
    drh = 0
    drcp = 0
    nur = 0
    for s,c in zip(subs,coeffs):
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
    ind = subs.index('O2') if 'O2' in subs else None
    coeff_O2 = coeffs[ind] if ind else 0

    ind = subs.index('H2O') if 'H2O' in subs else None
    coeff_H2O = coeffs[ind] if ind else 0

    ind = subs.index('CO2') if 'CO2' in subs else None
    coeff_CO2 = coeffs[ind] if ind else 0
    
    #Substances that will be approximated by Ntot
    inds = [subs.index(i) if i in subs else None for i in ['NO2','HNO3']]
    #There should be only one by equation...
    coeff_Ntot = sum([coeffs[ind] if ind else 0 for ind in inds])

    #Substances that will be approximated by Stot
    inds = [subs.index(i) if i in subs else None for i in ['COS','H2S','SO2','SO3','H2SO4']]
    #There should be only one by equation...
    coeff_Stot = sum([coeffs[ind] if ind else 0 for ind in inds])
    
    #P = {'S': 1, 'N': 1, 'C': 2000, 'T': 298.15}   <---------   Important!!!
    #The intercept of the equations is a function of the composition.
    intercept = lambda P, lgKp=_reactions[key]['lgK_p'], c_s=coeff_Stot, c_n=coeff_Ntot, c_c=coeff_CO2: -lgKp(P['T'])+c_s*np.log10(P['S'])+c_n*np.log10(P['N'])+c_c*np.log10(P['C'])
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
    
#Based on the specified rules, create a dictionary containing the lines specifying a graph
def _form_lines(Rules,_reactions=_reactions):
    Lines = dict()
    for key,rules in Rules.items():
        rules_f = list()
        for rule in rules:
            if rule[0]=='above':
                rules_f.append(lambda x, y, P, key=rule[1]: y > _reactions[key]['line']['y'](x,P))
            elif rule[0]=='below':
                rules_f.append(lambda x, y, P, key=rule[1]: y < _reactions[key]['line']['y'](x,P))
            elif rule[0]=='right of':
                rules_f.append(lambda x, y, P, key=rule[1]: x > _reactions[key]['line']['x'](y,P))
            elif rule[0]=='left of':
                rules_f.append(lambda x, y, P, key=rule[1]: x < _reactions[key]['line']['x'](y,P))
        
        Lines[key] = _reactions[key]['line'] | {'rules': rules, 'rules_f': rules_f}
        
    return Lines
    
#Create the lines dicts
_lines_Fe_O = _form_lines(_rules_Fe_O)
_lines_Fe_C = _form_lines(_rules_Fe_C)
_lines_Fe_N = _form_lines(_rules_Fe_N)
_lines_Fe_S = _form_lines(_rules_Fe_S)

#Add special rules
def _special_rule_Fe2O3_FeOOH(x,y,P):
    r1 = y < _reactions['Fe(NO3)2/NO2/FeO(OH)']['line']['y'](x,P)
    r2 = y < _reactions['Fe2O3/NO/Fe(NO3)2']['line']['y'](x,P)
    r3 = x > _math._intersection(_reactions['Fe(NO3)2/HNO3/FeO(OH)']['line'], _reactions['Fe(NO3)2/NO2/FeO(OH)']['line'], P)[0]
    
    return r1 or r2 or r3

def _special_rule_FeNO32_HNO3_FeOOH(x,y,P):
    i_1 = _math._intersection(_reactions['Fe(NO3)2/HNO3/FeO(OH)']['line'], _reactions['Fe2O3/NO2/Fe(NO3)2']['line'], P)
    i_2 = _math._intersection(_reactions['Fe2O3/HNO3/Fe(NO3)2']['line'], _reactions['Fe2O3/NO2/Fe(NO3)2']['line'], P)
    
    return i_1[0] < i_2[0] 
    
def _special_rule_Fe2O3_HNO3_FeNO32(x,y,P):
    i_1 = _math._intersection(_reactions['Fe(NO3)2/HNO3/FeO(OH)']['line'], _reactions['Fe2O3/NO2/Fe(NO3)2']['line'], P)
    i_2 = _math._intersection(_reactions['Fe2O3/HNO3/Fe(NO3)2']['line'], _reactions['Fe2O3/NO2/Fe(NO3)2']['line'], P)
    
    return i_1[0] >= i_2[0] 

_lines_Fe_N['Fe2O3/FeO(OH)']['rules_f'].append(_special_rule_Fe2O3_FeOOH)
_lines_Fe_N['Fe(NO3)2/HNO3/FeO(OH)']['rules_f'].append(_special_rule_FeNO32_HNO3_FeOOH)
_lines_Fe_N['Fe2O3/HNO3/Fe(NO3)2']['rules_f'].append(_special_rule_Fe2O3_HNO3_FeNO32)

def _special_rule_FeSO4H2O_FeSO4(x,y,P):  #('below', 'FeSO4/SO2/Fe2(SO4)3')
    r1 = y < _reactions['FeSO4/SO2/Fe2(SO4)3']['line']['y'](x,P)
    r2 = y < _reactions['FeSO4/S/Fe2(SO4)3']['line']['y'](x,P)
    
    return r1 or r2

_lines_Fe_S['FeSO4.H2O/FeSO4']['rules_f'].append(_special_rule_FeSO4H2O_FeSO4)

#Combine
_lines = {
    'O': _lines_Fe_O,
    'C': _lines_Fe_C,
    'N': _lines_Fe_N,
    'S': _lines_Fe_S,
}

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


