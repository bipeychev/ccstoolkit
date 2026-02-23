#!/usr/bin/python3

import numpy as np

#Get the angle between points p and q
def _angle(xy1: tuple, xy2: tuple):
    x1, y1 = xy1
    x2, y2 = xy2
    return np.atan2(y2-y1, x2-x1)

#Get the area of a polygon
def _polygon_area(face: list):
    area = 0.0
    for (x1, y1), (x2, y2) in zip(face, face[1:] + face[:1]):
        area += x1 * y2 - x2 * y1
    return area / 2

#Cross product
def _cross2d(x: tuple, y: tuple):
    return x[..., 0]*y[..., 1]-x[..., 1]*y[..., 0]

#Get the coordinates of the centroid
def _calculate_centroid(vertices: list):
    roll0 = np.roll(vertices, 0, axis=0)
    roll1 = np.roll(vertices, 1, axis=0)
    cross = _cross2d(roll0, roll1)
    area = 0.5 * np.sum(cross)
    return np.sum((roll0 + roll1) * cross[:, None], axis=0) / (6.0 * area)

#Format coordinates
def _format_xy(xy: tuple):
    x, y = xy
    return (round(float(x),6),round(float(y),6))     #Apparently lists are unhashable, but tuples are not
    
#Get the intersection between two lines
def _intersection(line_i: dict, line_j: dict, P: dict):  #a+b*x+c*y=0 
    a_i, b_i, c_i = line_i['coeffs'](P)
    a_j, b_j, c_j = line_j['coeffs'](P)
    
    det = b_i*c_j-b_j*c_i
    if abs(det) < 1e-9:               #parallel
        return None

    x = (c_i*a_j-c_j*a_i)/det
    y = -(b_i*a_j-b_j*a_i)/det
    return x, y
    
