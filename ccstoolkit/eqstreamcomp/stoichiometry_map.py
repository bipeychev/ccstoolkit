#!/usr/bin/python3

#===================================================================================================================
#---------------------------------------------------------------------------------------Import libraries
#===================================================================================================================
from . import _stoichiometry
import ccstoolkit.common._line_logic as _line_logic

_bounds = _stoichiometry.get_bounds()
_domain = _stoichiometry.get_domain()
_lines = _stoichiometry.get_lines()

#===================================================================================================================
#---------------------------------------------------------------------------------------Private methods
#===================================================================================================================
#Find the name of each region
def _parse_region_names(regions: dict):
    #Go through the regions
    for region in regions:
        #Remove the ids that are only digits. Those are the bounding box lines.
        bounding_lines = [i for i in region['bounds ids'] if not i.isdigit()]
        
        #Identify the unique substances on the bounding lines
        substances = list(set([i for line in bounding_lines for i in line.split('+')]))
        
        #In the instance where there is H2O, NO and NO2, the mixture is also at equilibrium with HNO2
        if {'H2O', 'NO', 'NO2'} <= set(substances):  #Is subset
            substances.append('HNO2')
        
        #Sort
        substances = sorted(substances)
        
        #Join into a label
        region['name'] = '+'.join(substances)
        
    return regions

#Get regions with names
def _get_regions_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
    return _parse_region_names(_line_logic._get_regions(lines, P, x_bounds, y_bounds))
    
#===================================================================================================================
#---------------------------------------------------------------------------------------Public methods
#===================================================================================================================
#A function that returns stoichiometric map for a given composition
def get_stoichiometry_map(P: dict):
    '''P = {
        'N/S': ratio of the nitrogen and sulphur concentrations
    }'''
    
    if not isinstance(P, dict):
        print('Wrong input!')
        return -1
        
        
    if {'N/S'} - set(P.keys()):
        print('Wrong input!')
        return -1
        
    for key, p in P.items():
        if not _domain[key]['min'] <= p <= _domain[key]['max']:
            print(f'Wrong input! {key} outside of range.')
            return -1
    
    regions = _get_regions_with_names(_lines, {'N': P['N/S'], 'S': 1}, _bounds['x'], _bounds['y'])
    
    #The region designated as 'COS+H2O+H2S+NO' is unexplored!
    return [region for region in regions if region['name']!='COS+H2O+H2S+NO']

