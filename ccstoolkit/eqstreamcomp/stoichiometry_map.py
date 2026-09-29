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
#Find the name of each face
def _parse_face_names(faces: dict):
	#Go through the faces
	for face in faces:
		#Remove the ids that are only digits. Those are the bounding box lines.
		bounding_lines = [i for i in face['bounds ids'] if not i.isdigit()]
		
		#Identify the unique substances on the bounding lines
		substances = list(set([i for line in bounding_lines for i in line.split('+')]))
		
		#In the instance where there is H2O, NO and NO2, the mixture is also at equilibrium with HNO2
		if {'H2O', 'NO', 'NO2'} <= set(substances):  #Is subset
			substances.append('HNO2')
		
		#Sort
		substances = sorted(substances)
		
		#Join into a label
		face['name'] = '+'.join(substances)
		
	return faces

#Get faces with names
def _get_faces_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
	return _parse_face_names(_line_logic._get_regions(lines, P, x_bounds, y_bounds))
	
#Find the name of each vertex
def _parse_vertex_names(vertices: dict):
	#Go through the vertices
	for vertex in vertices:
		#Find common substances in the ids		
		commons = set.intersection(*[set(i.split('+')) for i in vertex['ids']])
		
		#Join into a label
		#In the case of x=0, remove the hydrogen containing substances
		vertex['id'] = '+'.join(sorted(commons if vertex['xy'][0]!=0 else {i for i in commons if 'H' not in i}))
		
	return vertices
	
#Get the active vertices with names
def _get_vertices_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
	return _parse_vertex_names(_line_logic._get_vertices(lines, P, x_bounds, y_bounds))
	
#Verify input
def _verify_input(P: dict):
	'''P = {
		'N/S': ratio of the nitrogen and sulphur concentrations
	}'''
	
	if not isinstance(P, dict):
		print('Wrong input!')
		return False
		
		
	if {'N/S'} - set(P.keys()):
		print('Wrong input!')
		return False
		
	for key, p in P.items():
		if not _domain[key]['min'] <= p <= _domain[key]['max']:
			print(f'Wrong input! {key} outside of range.')
			return False
			
	return True
	
#===================================================================================================================
#---------------------------------------------------------------------------------------Public methods
#===================================================================================================================
#A function that returns the faces of the stoichiometric map for a given composition
def get_stoichiometry_map(P: dict):
	'''P = {
		'N/S': ratio of the nitrogen and sulphur concentrations
	}'''
	
	if _verify_input(P):
		regions = _get_faces_with_names(_lines, {'N': P['N/S'], 'S': 1}, _bounds['x'], _bounds['y'])
	
		#The region designated as 'COS+H2O+H2S+NO' is unexplored!
		return [region for region in regions if region['name']!='COS+H2O+H2S+NO']
	
#A function that returns the edges of the stoichiometric map for a given composition
def get_stoichiometry_map_lines(P: dict):
	'''P = {
		'N/S': ratio of the nitrogen and sulphur concentrations
	}'''
	
	if _verify_input(P):
		return _line_logic._get_edges(_lines, {'N': P['N/S'], 'S': 1}, _bounds['x'], _bounds['y'])
	
#A function that returns the vertices of the stoichiometric map for a given composition
def get_stoichiometry_map_points(P: dict):
	'''P = {
		'N/S': ratio of the nitrogen and sulphur concentrations
	}'''
	
	if _verify_input(P):
		return _get_vertices_with_names(_lines, {'N': P['N/S'], 'S': 1}, _bounds['x'], _bounds['y'])

