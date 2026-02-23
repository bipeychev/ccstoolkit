#!/usr/bin/python3

#===================================================================================================================
#---------------------------------------------------------------------------------------Import libraries
#===================================================================================================================
from . import _reactions
import ccstoolkit.common._line_logic as _line_logic

_bounds = _reactions.get_bounds()
_domain = _reactions.get_domain()
_lines = _reactions.get_lines()

#===================================================================================================================
#---------------------------------------------------------------------------------------Private methods
#===================================================================================================================
#Find the name of each region
def _parse_region_names(regions: dict):
    #Go through the regions
    for region in regions:
        #Remove the ids that are only digits. Those are the bounding box lines.
        bounding_lines = [i for i in region['bounds ids'] if not i.isdigit()]
        
        #Split the id, e.g. 'H2S/S+NO' -> {'H2S','S','NO'}
        j = [set([k for j in i.split('/') for k in j.split('+')]) for i in bounding_lines]
        
        #The name of the region are the common substances for all lines, e.g. {'S','NO'}
        region['name'] = j[0].intersection(*j[1:])
        
        #Combine
        region['name'] = '+'.join(sorted(region['name']))
        
    return regions

#Get regions with names
def _get_regions_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
    return _parse_region_names(_line_logic._get_regions(lines, P, x_bounds, y_bounds))
    
#===================================================================================================================
#---------------------------------------------------------------------------------------Public methods
#===================================================================================================================
#A function that returns the stability map for a given composition
def get_stability_map(P: dict):
    '''P = {
        'S': total sulphur concentration in [mM], 
        'N': total nitrogen concentration in [mM], 
        'CO2': activity of CO2 in [mM], 
        'T': temperature in [K]
    }'''
    
    if not isinstance(P, dict):
        print('Wrong input!')
        return -1
        
    #Add default values, if not specified
    for key, default in [('CO2', 2e3),('T', 298.15)]:
        if not key in P:
            P = P | {key: default}
        
    if {'S','N'} - set(P.keys()):
        print('Wrong input!')
        return -1
        
    for key, p in P.items():
        if not _domain[key]['min'] <= p <= _domain[key]['max']:
            print(f'Wrong input! {key} outside of range.')
            return -1
    
    return _get_regions_with_names(_lines, P, _bounds['x'], _bounds['y'])

