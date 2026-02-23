#!/usr/bin/python3
    
#Properties of substances
_substances_TD_data = {
    #substance: {dfg in kJ/mol, dfh in kJ/mol, cp in J/mol/K, solid ? True : False}
    'O2': {'dfg': 0, 'dfh': 0, 'cp': 29.4, 'solid': False},
    'S': {'dfg': 0, 'dfh': 0, 'cp': 22.6,  'solid': True},
    'Fe': {'dfg': 0, 'dfh': 0, 'cp': 25.1, 'solid': True},
    #--------------------------------------
    'H2O': {'dfg': -228.6, 'dfh': -241.8, 'cp': 33.6,  'solid': False},
    'CO2': {'dfg': -394.4, 'dfh': -393.5, 'cp': 37.1,  'solid': False},
    'CO': {'dfg': -137.2, 'dfh': -110.5, 'cp': 29.1,  'solid': False},
    'NO': {'dfg': 87.6, 'dfh': 91.3, 'cp': 29.9,  'solid': False},
    'NO2': {'dfg': 51.3, 'dfh': 33.2, 'cp': 37.2,  'solid': False},
    'H2S': {'dfg': -33.4, 'dfh': -20.6, 'cp': 34.2,  'solid': False},
    'COS': {'dfg': -169.2, 'dfh': -142.0, 'cp': 41.5,  'solid': False},
    'SO2': {'dfg': -300.1, 'dfh': -296.8, 'cp': 39.9,  'solid': False},
    'SO3': {'dfg': -371.1, 'dfh': -395.7, 'cp': 50.7,  'solid': False},
    'H2SO4': {'dfg': -653.31, 'dfh': -735.13, 'cp': 83.7,  'solid': False},
    'HNO2': {'dfg': -46.0, 'dfh': -79.5, 'cp': 45.6,  'solid': False},
    'HNO3': {'dfg': -73.5, 'dfh': -133.9, 'cp': 54.1,  'solid': False},
    #--------------------------------------
    'FeO': {'dfg': -251.46, 'dfh': -271.96, 'cp': 49.915, 'solid': True},
    'Fe3O4': {'dfg': -1015.46, 'dfh': -1118.4, 'cp': 143.4, 'solid': True},
    'Fe2O3': {'dfg': -742.24, 'dfh': -824.25, 'cp': 103.9, 'solid': True},
    'FeCO3': {'dfg': -666.72, 'dfh': -740.57, 'cp': 82.1, 'solid': True},
    'Fe(OH)3': {'dfg': -696.5, 'dfh': -823.0, 'cp': 101.671, 'solid': True},
    'Fe(OH)2': {'dfg': -486.5, 'dfh': -569.0, 'cp': 97.069, 'solid': True},
    'FeO(OH)': {'dfg': -489.8, 'dfh': -561.9, 'cp': 74.5, 'solid': True},
    'FeS2': {'dfg': -166.94, 'dfh': -178.24, 'cp': 62.2, 'solid': True},
    'FeS': {'dfg': -100.42, 'dfh': -99.998, 'cp': 50.5, 'solid': True},
    'FeSO4': {'dfg': -825.08, 'dfh': -928.43, 'cp': 100.6, 'solid': True},
    'Fe2(SO4)3': {'dfg': -2263.1, 'dfh': -2581.5, 'cp': 264.722, 'solid': True},
    'FeSO4.H2O': {'dfg': -1081.2, 'dfh': -1243.69, 'cp': 0, 'solid': True},                 #Can't find data for cp
    'FeSO4.7H2O': {'dfg': -2510.274, 'dfh': -3014.572, 'cp': 394.5, 'solid': True},
    'Fe(NO3)2': {'dfg': -288, 'dfh': -497.9, 'cp': 150, 'solid': True}
}

def get_substances_TD_data():
    return _substances_TD_data
    
