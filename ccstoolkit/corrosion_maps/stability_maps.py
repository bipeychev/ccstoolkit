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
#Find the name of each face
def _parse_face_names(faces: dict):
	#Go through the faces
	for face in faces:
		#Remove the ids that are only digits. Those are the bounding box lines.
		bounding_lines = [i for i in face['bounds ids'] if not i.isdigit()]
		
		#Split the id, e.g. 'FeS/H2S/FeS2' -> {'FeS','H2S','FeS2'}
		j = [set(i.split('/')) for i in bounding_lines]
		
		#The name of the face are the common substances for all lines, e.g. {'FeS'}
		face['name'] = j[0].intersection(*j[1:])
		
	#Edge faces might not have enough lines to single out one single common substance.
	all_names = [face['name'] for face in faces]
	#Go through all faces
	for face in faces:
		this_name = face['name']
		other_names = [name for name in all_names if name!=this_name]
		
		#If more than one elements are present in 'name', i.e. the algorithm wasn't able to identify a face
		if len(this_name)>1:
			#Remove common (common with other faces) names and obvious non-corrosion products
			for name in other_names+[{'CO2'},{'HNO3'},{'H2S'}]:
				this_name -= name   
				
	#Just in case there is a face left without a name, ignore it (it's called problem solving!)
	faces = [face for face in faces if face['name']]
	
	#Flatten the names
	for face in faces:
		face['name'] = list(face['name'])[0]
		
	return faces

#Get faces with names
def _get_faces_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
	return _parse_face_names(_line_logic._get_regions(lines, P, x_bounds, y_bounds))
	
#Find the name of each edge
def _parse_edge_names(edges: dict):
	#Go through the edges
	for edge in edges:
		#Format the id	
		#There are non-corrosion product substances for uniqueness
		edge['id'] = ' = '.join(sorted([i for i in edge['id'].split('/') if 'Fe' in i]))
		
	return edges

#Get the active edges with names
def _get_edges_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
	return _parse_edge_names(_line_logic._get_edges(lines, P, _bounds['x'], _bounds['y']))
	
#Find the name of each vertex
def _parse_vertex_names(vertices: dict):
	#Go through the vertices
	for vertex in vertices:

		vertex['id'] = ' = '.join(sorted({j for i in vertex['ids'] for j in i.split('/') if 'Fe' in j}))
		
	return vertices
	
#Get the active vertices with names
def _get_vertices_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
	return _parse_vertex_names(_line_logic._get_vertices(lines, P, x_bounds, y_bounds))
	
#Add default values
def _add_default_values(P: dict):
	for key, default in [('CO2', 2e3),('T', 298.15)]:
		if key not in P:
			P = P | {key: default}
			
	return P

#Verify input
def _verify_input(P: dict):
	'''P = {
		'S': total sulphur concentration in [mM], 
		'N': total nitrogen concentration in [mM], 
		'CO2': activity of CO2 in [mM], 
		'T': temperature in [K]
	}'''
	
	if not isinstance(P, dict):
		print('Wrong input!')
		return False
		
	if {'S','N'} - set(P.keys()):
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
#A function that returns all graphs for a given composition
def get_stability_maps(P: dict):
	'''P = {
		'S': total sulphur concentration in [mM], 
		'N': total nitrogen concentration in [mM], 
		'CO2': activity of CO2 in [mM], 
		'T': temperature in [K]
	}'''
	
	#Add default values, if not specified
	P = _add_default_values(P)
	
	#Get regions	
	if _verify_input(P):
		regions = {key: _get_faces_with_names(lines, P, _bounds['x'], _bounds['y']) for key, lines in _lines.items()}
		
		return regions
		
#A function that returns the edges of the stability map for a given composition
def get_stability_maps_lines(P: dict):
	'''P = {
		'S': total sulphur concentration in [mM], 
		'N': total nitrogen concentration in [mM], 
		'CO2': activity of CO2 in [mM], 
		'T': temperature in [K]
	}'''
	
	#Add default values, if not specified
	P = _add_default_values(P)
	
	#Get regions	
	if _verify_input(P):
		return {key: _get_edges_with_names(lines, P, _bounds['x'], _bounds['y']) for key, lines in _lines.items()}

#A function that returns the vertices of the stability map for a given composition
def get_stability_maps_points(P: dict):
	'''P = {
		'S': total sulphur concentration in [mM], 
		'N': total nitrogen concentration in [mM], 
		'CO2': activity of CO2 in [mM], 
		'T': temperature in [K]
	}'''
	
	#Add default values, if not specified
	P = _add_default_values(P)
	
	#Get regions	
	if _verify_input(P):
		return {key: _get_vertices_with_names(lines, P, _bounds['x'], _bounds['y']) for key, lines in _lines.items()}

