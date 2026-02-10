#!/usr/bin/python3

#===================================================================================================================
#---------------------------------------------------------------------------------------Import libraries
#===================================================================================================================
import numpy as np

from . import _math
from . import _reactions

_bounds = _reactions.get_bounds()
_lines = _reactions.get_lines()

#===================================================================================================================
#---------------------------------------------------------------------------------------Private methods
#===================================================================================================================
    
#Get all intersections between a list of lines
def _all_intersections(lines: dict, P: dict):
    ids = list(lines.keys())
    points = []

    for i in range(len(ids)):
        for j in range(i+1, len(ids)):
            p = _math._intersection(lines[ids[i]], lines[ids[j]], P)
            if p is not None:
                points.append(p)
    return points

#Get only the active part of a line constrained by x_bounds, y_bounds and applying the line['rules']
def _clip_line(line: dict, intersections: list, P: dict, x_bounds: tuple, y_bounds: tuple):
    x_min, x_max = x_bounds
    y_min, y_max = y_bounds
    a, b, c = line['coeffs'](P)
    
    #Collect coordinates where this line intersects the others
    t_vals = []                              #t can be x or y depending on the orientation
    for p in intersections:
        if p is None:
            continue
        x, y = p

        if abs(a+b*x+c*y)<1e-9:               #Satisfies the equation = lies on the line
            t_vals.append(y if line['vertical'] else x)
            
    #For angled lines, add intersections with the bounding box
    if not line['vertical'] and not line['horizontal']:
        t_vals.append(line['y'](x_min, P))
        t_vals.append(line['y'](x_max, P))
        t_vals.append(line['x'](y_min, P))
        t_vals.append(line['x'](y_max, P))
            
    #Add bounding box limits
    t_vals += (y_bounds if line['vertical'] else x_bounds)
    
    #Remove intersection points ouside the bounds
    if line['vertical']:
        t_vals = [t for t in t_vals if y_min <= round(t,6) <= y_max and x_min <= round(line['x'](t, P),6) <= x_max]
    else:
        t_vals = [t for t in t_vals if x_min <= round(t,6) <= x_max and y_min <= round(line['y'](t, P),6) <= y_max]
    
    #Sort
    t_vals = sorted(t_vals)
    
    #Generate segments between consecutive points
    segments = list()     
    for t0, t1 in zip(t_vals[:-1], t_vals[1:]):
        tmid = 0.5*(t0+t1)

        if line['vertical']:
            xmid = line['x'](tmid, P)
            ymid = tmid
        else:
            xmid = tmid
            ymid = line['y'](tmid, P)
            
        #If all rules are satisfied, save this segment
        if np.all([rule(xmid,ymid,P) for rule in line['rules_f']]):
            if line['vertical']:
                segments.append(((xmid, t0), (xmid, t1)))
            else:
                segments.append((
                    (t0, line['y'](t0, P)),
                    (t1, line['y'](t1, P))
                ))

    return segments
    
#Get all active segments
def _get_active_lines(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
    #Get all intersection points
    intersections = _all_intersections(lines, P)

    active_lines = list()
    #For each line
    for key, line in lines.items():
        #Check if active
        if 'active' in line:
            if not line['active'](P):
                continue
        
        #Get the active segments of the line
        segments = _clip_line(line, intersections, P, x_bounds, y_bounds)
        #The segments are sorted by x or y!
        
        #If there are active segments
        if segments:
            #Combine the segment into one.
            #Works because: all lines are straight and the segments are sorted
            #The points need to be rounded, otherwise truncation errors make the match problematic!
            active_line = dict()
            active_line['id'] = key                               #The id of the line
            active_line['p0'] = _math._format_xy(segments[0][0])         #One end point
            active_line['p1'] = _math._format_xy(segments[-1][1])        #The other end point
            active_lines.append(active_line)
            
    return active_lines

#Get the regions defined by a list of lines and x,y bounds
def _get_regions(lines: dict, P: dict, x_bounds: tuple, y_bounds: tuple):
    #---------------------------Add bounding box segments---------------------------
    #The corners of the bounding box
    corners = [(x, y) for x in x_bounds for y in y_bounds]

    #The active lines
    active_lines = _get_active_lines(lines, P, x_bounds, y_bounds)

    #Get all intersection points
    active_intersections = [p for line in active_lines for p in [line['p0'],line['p1']]]
    active_intersections = set(active_intersections)

    #Get the intersection points that lie on the bounding box
    intersections_box = [(x,y) for x,y in active_intersections if x in x_bounds or y in y_bounds]

    #Separate the edges of the bounding box    #It's important that they are sorted!
    left_right = [sorted([(x,y) for x,y in intersections_box+corners if x==xlim],key=lambda p: p[1]) for xlim in x_bounds]
    top_bottom = [sorted([(x,y) for x,y in intersections_box+corners if y==ylim],key=lambda p: p[0]) for ylim in y_bounds]
    box_edges = left_right+top_bottom

    #Form line segments from box_edges
    lines_box = list()
    for cnt, edge in enumerate(box_edges):
        for cnt1, p in enumerate(edge[:-1]):
            line = {'id': str(cnt)+str(cnt1), 'p0': p, 'p1': edge[cnt1+1]}
            lines_box.append(line)

    #---------------------------Find the faces---------------------------
    #Walk the edges starting from a point and always turning in the same direction.
    #When the point repeats, a face has been enclosed.
    
    #Define inner and outer edges for each line segment (inner/outer = top/bottom = left/right)
    edges = list()
    for line in active_lines+lines_box:
        edges.append({'id': line['id'], 'from': line['p0'], 'to': line['p1'], 'used': False})
        edges.append({'id': line['id'], 'from': line['p1'], 'to': line['p0'], 'used': False})

    #Get all vertices
    vertices = dict()
    for edge in edges:
        vertice = edge['from']

        if vertice not in vertices:
            vertices[vertice] = []

        vertices[vertice].append(edge)

    #For each vertice, sort the outgoing edges by angle.
    #This ensures the algoritm transverses the graph always turning in the same direction (left).
    for vertice, out_edges in vertices.items():
        out_edges.sort(key=lambda e: _math._angle(vertice, e['to']))
        
    #Get all regions
    regions = list()
    #Starting from an edge
    for edge in edges:
        #Only if the edge is active
        if edge['used']:
            continue

        #face = list()
        face = {'ids': [], 'points': []}
        current = edge
        
        #Add all 'from' points untill the loop closes
        while True:
            current['used'] = True                             #Remove this edge from rotation
            
            face['ids'].append(current['id'])
            face['points'].append(current['from'])             #Starting point

            towards = current['to']                            #End point
            outs = vertices[towards]                           

            #Find the index of the reverse edge current['to']->current['from']
            for i, oe in enumerate(outs):
                if oe['to'] == current['from']:
                    idx = i
                    break

            #Turn left: previous edge in CCW order
            current = outs[(idx - 1) % len(outs)]

            #Check if loop is closed
            if current is edge:
                break

        regions.append(face)

    #---------------------------Format the dictionary---------------------------
    #Close the loops for plotting
    regions = [{'ids': region['ids'], 'points': region['points']+[region['points'][0]]} for region in regions]
    
    #Calculate the area of each face
    for region in regions:
        region['area'] = _math._polygon_area(region['points'])
        
    #One gets the inner faces + the outer (infinite) face of the bounding box.
    #The latter will always have the largest area + a different sign of the area.
    inner_regions = [region for region in regions if region['area']>0]
    
    #Determine the name of each region
    for region in inner_regions:
        #Remove the ids that are only digits. Those are the bounding box lines.
        ids = [i for i in region['ids'] if not i.isdigit()]
        #Split the id, e.g. 'FeS/H2S/FeS2' -> {'FeS','H2S','FeS2'}
        j = [set(i.split('/')) for i in ids]
        #The name of the region is the common element for all lines
        region['name'] = j[0].intersection(*j[1:])
    
    #Edge regions might not have enough lines to single out one single common substance.
    all_names = [region['name'] for region in inner_regions]
    #Go through all regions
    for region in inner_regions:
        this_name = region['name']
        other_names = [name for name in all_names if name!=this_name]
        #If more than one elements are present in 'name', i.e. the algorithm wasn't able to identify a region
        if len(this_name)>1:
            #Remove common (common with other regions) names and obvious non-names
            for name in other_names+[{'CO2'},{'HNO3'},{'H2S'}]:
                this_name -= name
    
    #Just in case there is a region without a name, remove it
    inner_regions = [region for region in inner_regions if region['name']]
    
    #Flatten the names
    for region in inner_regions:
        region['name'] = list(region['name'])[0]

    #Calculate centroids
    for region in inner_regions:
        region['centroid'] = _math._calculate_centroid(region['points'])
        
    #Sort the regions and keys
    key_order = ['name','area','centroid','points']
    inner_regions = sorted([{key: region[key] for key in key_order if key in region} for region in inner_regions], key=lambda d: d["name"])
        
    return inner_regions

#===================================================================================================================
#---------------------------------------------------------------------------------------Public methods
#===================================================================================================================
#A function that returns all graphs for a given composition
def get_maps(P: dict):
    '''P = {
        'S': total sulphur concentration in [mM], 
        'N': total nitrogen concentration in [mM], 
        'C': activity of CO2 in [mM], 
        'T': temperature in [K]
    }'''
    
    Regions = {key: _get_regions(lines, P, _bounds['x'], _bounds['y']) for key, lines in _lines.items()}
    
    return Regions


