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
        
        #Split the id, e.g. 'FeS/H2S/FeS2' -> {'FeS','H2S','FeS2'}
        j = [set(i.split('/')) for i in bounding_lines]
        
        #The name of the region are the common substances for all lines, e.g. {'FeS'}
        region['name'] = j[0].intersection(*j[1:])
        
    #Edge regions might not have enough lines to single out one single common substance.
    all_names = [region['name'] for region in regions]
    #Go through all regions
    for region in regions:
        this_name = region['name']
        other_names = [name for name in all_names if name!=this_name]
        
        #If more than one elements are present in 'name', i.e. the algorithm wasn't able to identify a region
        if len(this_name)>1:
            #Remove common (common with other regions) names and obvious non-corrosion products
            for name in other_names+[{'CO2'},{'HNO3'},{'H2S'}]:
                this_name -= name   
                
    #Just in case there is a region left without a name, ignore it
    regions = [region for region in regions if region['name']]
    
    #Flatten the names
    for region in regions:
        region['name'] = list(region['name'])[0]
        
    return regions

#Get regions with names
def _get_regions_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
    return _parse_region_names(_line_logic._get_regions(lines, P, x_bounds, y_bounds))
    
#===================================================================================================================
#---------------------------------------------------------------------------------------Public methods
#===================================================================================================================
#A function that returns all graphs for a given composition
def get_stability_maps(P: dict):
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
        if key not in P:
            P = P | {key: default}
    
    if {'S','N'} - set(P.keys()):
        print('Wrong input! Insufficient number of concentrations provided.')
        return -1
        
    for key, p in P.items():
        if not _domain[key]['min'] <= p <= _domain[key]['max']:
            print(f'Wrong input! {key} outside of range.')
            return -1
    
    regions = {key: _get_regions_with_names(lines, P, _bounds['x'], _bounds['y']) for key, lines in _lines.items()}
    
    return regions


