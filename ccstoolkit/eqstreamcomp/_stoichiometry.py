#!/usr/bin/python3

#===================================================================================================================
#---------------------------------------------------------------------------------------Import
#===================================================================================================================
import re
import numpy as np
import numpy.linalg as linalg

import ccstoolkit.common._line_logic as _line_logic
#===================================================================================================================
#---------------------------------------------------------------------------------------Data
#===================================================================================================================
_bounds = {
    'x': (0, 9),
    'y': (0, 9)
}

_domain = {
    'S': {'min': 0.015, 'max': 4},                                                       #[mol/m^3]
    'N': {'min': 0.015, 'max': 4},                                                       #[mol/m^3]
    'N/S': {'min': 0.050, 'max': 4}                                                      #[-]
}

#The stoichiometric lines that are going to be plot
_unformatted_lines = {
    'COS+H2S+NO': dict(),
    'H2O+H2S+NO': dict(),
    'H2O+H2SO4+HNO3': dict(),
    'H2O+H2SO4+NO': dict(),
    'H2O+H2SO4+NO2': dict(),
    'H2O+NO+S': dict(),
    'H2O+NO+SO2': dict(),
    'H2S+NO+S': dict(),
    'H2SO4+HNO3+NO2': dict(),
    'H2SO4+HNO3+O2': dict(),
    'H2SO4+NO+NO2': dict(),
    'H2SO4+NO+SO2': dict(),
    'H2SO4+NO+SO3': dict(),
    'H2SO4+NO2+O2': dict(),
    'H2SO4+NO2+SO3': dict()
}

#Specify each graph as a set of rules bounding the domain of each line
#Only lines with a rule will be included in the graph (the rules can be an empty list)
_rules = {
    'COS+H2S+NO': [('below','H2S+NO+S')],
    'H2O+H2S+NO': [('above','H2S+NO+S')],
    'H2O+H2SO4+HNO3': [('right of', 'H2SO4+HNO3+O2')],
    'H2O+H2SO4+NO': [('right of', 'H2SO4+NO+NO2')],
    'H2O+H2SO4+NO2': [('right of', 'H2SO4+NO2+O2')],
    'H2O+NO+S': [],
    'H2O+NO+SO2': [],
    'H2S+NO+S': [('left of','H2O+H2S+NO')],
    'H2SO4+HNO3+NO2': [('above', 'H2O+H2SO4+NO2'), ('below', 'H2O+H2SO4+HNO3')],
    'H2SO4+HNO3+O2': [('above', 'H2O+H2SO4+HNO3')],
    'H2SO4+NO+NO2': [('above', 'H2SO4+NO+SO3'),('below','H2SO4+NO2+SO3')],
    'H2SO4+NO+SO2': [('below', 'H2SO4+NO+SO3')],
    'H2SO4+NO+SO3': [('left of', 'H2SO4+NO+NO2')],
    'H2SO4+NO2+O2': [('above', 'H2O+H2SO4+NO2')],
    'H2SO4+NO2+SO3': [('left of', 'H2SO4+NO2+O2')]
}

#===================================================================================================================
#---------------------------------------------------------------------------------------Calculate/format/etc.
#===================================================================================================================
cti = lambda x: int(x) if x else 1

#Format the lines
for line, i in _unformatted_lines.items():
    #The substances on each line
    substances = line.split('+')
    
    #The stoichiometric coefficients
    H = np.array([cti(re.findall(r'H(\d*)', s)[0]) if 'H' in s else 0 for s in substances])
    N = np.array([cti(re.findall(r'N(\d*)', s)[0]) if 'N' in s else 0 for s in substances])
    O = np.array([(-1 if s=='COS' else 1)*cti(re.findall(r'O(\d*)', s)[0]) if 'O' in s else 0 for s in substances])
    S = np.array([cti(re.findall(r'S(\d*)', s)[0]) if 'S' in s else 0 for s in substances])
    
    i['H'] = H
    i['N'] = N
    i['O'] = O
    i['S'] = S
    
    #Define the line equation
    #P = {'N': 1, 'S': 1}   <---------   Important!!!
    if linalg.det(np.array([H,N,S]))!=0:
        matrix = np.array([H,N,S])
        vector = O @ linalg.inv(matrix)      #Matrix multiply
        
        #The intercept of the equations is a function of the composition.
        intercept = lambda P, vector=vector: vector[1:] @ np.array([P['N'],P['S']])
        slope = vector[0]
        
        #i['O @ coeff_matrix'] = vector
        i['equation'] = 'c[0](P)+c[1]*c_H+c[2]*c_O=0'
        i['coeffs'] = lambda P, c=(intercept,slope): (c[0](P),c[1],-1)
        
        #---------------------------Equations in the form c_O=f(c_H) and c_H=f(c_O)---------------------------
        i['x'] = lambda y, P, c=(intercept,slope): (y-c[0](P))/c[1]
        i['y'] = lambda x, P, vector=vector: vector @ np.array([x,P['N'],P['S']])
        
        #---------------------------Is the line horizontal/vertical---------------------------
        i['vertical'] = False
        i['horizontal'] = slope == 0
    else:                                    #When the determinant is 0, the line is vertical
        matrix = np.array([O,N,S])
        vector = H @ linalg.inv(matrix)      #Matrix multiply
        
        #The intercept of the equations is a function of the composition.
        intercept = lambda P, vector=vector: vector[1:] @ np.array([P['N'],P['S']])
        slope = vector[0]
        
        i['equation'] = 'c[0](P)+c[1]*c_H+c[2]*c_O=0'
        i['coeffs'] = lambda P, c=(intercept,slope): (c[0](P),-1,c[1])
        
        #---------------------------Equations in the form c_O=f(c_H) and c_H=f(c_O)---------------------------
        i['x'] = lambda y, P, vector=vector: vector @ np.array([y,P['N'],P['S']])
        i['y'] = lambda x, P, c=(intercept,slope): (x-c[0](P))/c[1]
        
        #---------------------------Is the line horizontal/vertical---------------------------
        i['vertical'] = True
        i['horizontal'] = False
        
#Create the lines dicts
_lines = _line_logic._form_lines(_rules,_unformatted_lines)

#===================================================================================================================
#---------------------------------------------------------------------------------------Export
#===================================================================================================================
def get_domain():
	return _domain

def get_bounds():
    return _bounds
    
def get_lines():
    return _lines
