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
		
		#Split the id, e.g. 'H2S/S+NO' -> {'H2S','S','NO'}
		j = [set([k for j in i.split('/') for k in j.split('+')]) for i in bounding_lines]
		
		#The name of the face are the common substances for all lines, e.g. {'S','NO'}
		face['name'] = j[0].intersection(*j[1:])
		
		#Combine
		face['name'] = '+'.join(sorted(face['name']))
		
	return faces

#Get faces with names
def _get_faces_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
	return _parse_face_names(_line_logic._get_regions(lines, P, x_bounds, y_bounds))
	
#Find the name of each edge
def _parse_edge_names(edges: dict):
	#Go through the edges
	for edge in edges:
		#Format the id
		edge['id'] = ' = '.join(sorted(edge['id'].split('+')[0].split('/')))			#The equistable tuple is infront by design
		
	return edges

#Get the active edges with names
def _get_edges_with_names(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
	return _parse_edge_names(_line_logic._get_edges(lines, P, _bounds['x'], _bounds['y']))
	
#Find the name of each vertex
def _parse_vertex_names(vertices: dict):
	#Go through the vertices
	for vertex in vertices:
		#The ids of the edges
		ids = [i.split('+')[0] for i in vertex['ids']]
		
		#If only one, this is an intersection with the bounding box
		if len(ids)==1:
		
			vertex['id'] = ids[0].replace('/',' = ')
			
		#If multiple
		else:
		
			#Create a container
			groups = list()
			
			#For each air of equistable substances
			for pair in [i.split('/') for i in ids]:
				#For each group
				for group in groups:
					#If either substance is already in a group, add the pair to the group
					if any(x in pair for x in group):
						[group.append(x) for x in pair]
						break
						
				#If there was nothing common found
				else:
					#Create a new group
					groups.append([*pair])
					
			vertex['id'] = '\n'.join(sorted([' = '.join(sorted(set(group))) for group in groups]))
		
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
#A function that returns the faces of the stability map for a given composition
def get_stability_map(P: dict):
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
		return _get_faces_with_names(_lines, P, _bounds['x'], _bounds['y'])
		
#A function that returns the edges of the stability map for a given composition
def get_stability_map_lines(P: dict):
	'''P = {
		'S': total sulphur concentration in [mM], 
		'N': total nitrogen concentration in [mM], 
		'CO2': activity of CO2 in [mM], 
		'T': temperature in [K]
	}'''
	
	#Add default values, if not specified
	P = _add_default_values(P)
	
	if _verify_input(P):
		return _get_edges_with_names(_lines, P, _bounds['x'], _bounds['y'])
		
#A function that returns the vertices of the stability map for a given composition
def get_stability_map_points(P: dict):
	'''P = {
		'S': total sulphur concentration in [mM], 
		'N': total nitrogen concentration in [mM], 
		'CO2': activity of CO2 in [mM], 
		'T': temperature in [K]
	}'''
	
	#Add default values, if not specified
	P = _add_default_values(P)
	
	if _verify_input(P):
		return _get_vertices_with_names(_lines, P, _bounds['x'], _bounds['y'])

