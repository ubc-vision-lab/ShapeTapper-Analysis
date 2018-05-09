#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 17:02:22 2018

Generates "heat map" plot of touch points data for each shape

Heat maps are output to the "heat_maps/" subdirectory

@author: Jamie Dunkle
"""

import cv2
import os, errno
import numpy as np
from ShapeObj import Shape
from ShapeIO import ShapeIO
from map_shape_pairs import rotate_points


################# Function Definitions #########################################################################
def plotHeatMap(shape, out_path, patient, cond=None) :

    print "Plotting {0}".format(shape.name)

    if patient is None :
        print "Error in plotHeatMap() : no patient name specified."
        return

    if shape.observed is None :
        print "Error in plotHeatMap() : no observed touch points found for {0}. Please check data files.".format(shape.name)
        return

    out_path = os.path.join(out_path, patient, "figures", "heat_maps")
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    observed = shape.observed
    medial_axis = shape.medial_axis.astype(np.uint32)
    edge_points = shape.edge_points.astype(np.uint32)
    centroid = shape.centroid.astype(np.uint32)

    img_bounds = np.array( [ [0,0], [shape.dims[0], 0] , [shape.dims[0],shape.dims[1]], [0,shape.dims[1]] ])
    circ_cent, circ_rad = cv2.minEnclosingCircle(img_bounds)
    circ_bd_x = [np.floor(circ_cent[1] - circ_rad), np.ceil(circ_cent[1] + circ_rad)]
    circ_bd_y = [np.floor(circ_cent[0] - circ_rad), np.ceil(circ_cent[0] + circ_rad)]
    canvas_size = circ_bd_x[1] - circ_bd_x[0]
    canvas_offset = np.array([int(circ_bd_y[0]),int(circ_bd_x[0])]) #(x,y)
    
    canvas = np.zeros((int(canvas_size),int(canvas_size)), dtype=np.uint8)
    heat_map = np.zeros((int(canvas_size),int(canvas_size),4), dtype=np.uint8)
    
    circ_cent = np.array([int(circ_cent[0]),int(circ_cent[1])])
    cv2.circle( canvas, tuple(circ_cent-canvas_offset), int(circ_rad), 0, thickness=-1 )
    
    for op in observed :
        o = np.array([int(op[0]),int(op[1])]) - np.flip(canvas_offset,axis=0) #(y,x)
        overlay = canvas.copy()
        rad = (int) (np.min(shape.dims[0:2]) / 10) # point radius
        for i in range(1, rad+1) :
            alpha = np.true_divide(1,3*rad)
            cv2.circle( overlay, tuple(o), i, 255, thickness=-1, lineType=cv2.LINE_AA )
            cv2.addWeighted(overlay, alpha, canvas, 1 - alpha, 0, canvas)  

    heat_map = cv2.applyColorMap(canvas, cv2.COLORMAP_JET)

    b_channel, g_channel, r_channel = cv2.split(heat_map)
    alpha_channel = np.ones(b_channel.shape, dtype=b_channel.dtype) * 255 # creating a dummy alpha channel image.
    heat_map = cv2.merge((b_channel, g_channel, r_channel, alpha_channel))
    
    heat_map[heat_map[:,:,0]==128] = [0,0,0,0] # get blue pixels (bottom of jet color map) & convert to transparent

    img_mask = np.zeros((int(canvas_size),int(canvas_size),4), dtype=np.uint8)
    cv2.fillConvexPoly(img_mask, edge_points - np.flip(canvas_offset,axis=0) , (255,255,255,255))
    img_mask = cv2.bitwise_or(img_mask,heat_map)
    heat_map[heat_map[:,:,3]==0] = [255,255,255,255]
    heat_map = cv2.bitwise_and(heat_map,img_mask)

    line_width = np.max(shape.dims[0:2])/95 # set line width according to image dimensions

    for mp in medial_axis :
        m = mp - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( heat_map, tuple(m), line_width, (255,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
    
    c = centroid[0] - np.flip(canvas_offset,axis=0) #(y,x)
    cv2.circle( heat_map, tuple(c), line_width*3, (0,255,0,255), thickness=-1, lineType=cv2.LINE_AA )  

    # heat_map[heat_map[:,:,3]==0] = [0,0,0,255] # convert shape interior to transparent

    for ep in edge_points :
        e = ep - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( heat_map, tuple(e), line_width, (0,0,0,255), thickness=-1, lineType=cv2.LINE_AA )

    if shape.pair_mapping is not None :
        out_fname = "_".join((shape.pair_mapping, "to", shape.name, "Patient", patient, 'heat_map.png'))
    else :
        out_fname = "_".join((shape.name, "Patient", patient, 'heat_map.png'))

    cv2.imwrite( os.path.join(out_path,out_fname), heat_map )



def plotMedialAxis(shape, out_path, patient=None, cond=None) :

    print "Plotting {0}".format(shape.name)

    out_path = os.path.join(out_path, "Shapes", "figures", "medial_axis")
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    img = shape.img.copy()
    medial_axis = shape.medial_axis.astype(np.uint32)
    edge_points = shape.edge_points.astype(np.uint32)
    centroid = shape.centroid.astype(np.uint32)
    
    line_width = np.max(img.shape[0:2])/190 # set line width according to image dimensions
    for mp in medial_axis :
       cv2.circle( img, tuple(mp), line_width, (255,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
   
    cv2.circle( img, tuple(centroid[0]), line_width*4, (0,255,0,255), thickness=-1, lineType=cv2.LINE_AA )  
    
    for ep in edge_points :
        cv2.circle( img, tuple(ep), line_width, (0,0,0,255), thickness=-1, lineType=cv2.LINE_AA )

    # heat_map[heat_map[:,:,3]==0] = [0,0,0,255] # convert shape interior to transparent

    out_fname = "_".join((shape.name, '_ma.png'))
    cv2.imwrite( os.path.join(out_path, out_fname), img )



def plotMappedShapes(shape_to, out_path, patient, cond=None) :

    if shape_to.pair_mapping is None :
        print "Error in plotMappedShapes() : list of single shape names. Please use a list of shape pairs."
        return

    if patient is None :
        print "Error in plotMappedShapes() : no patient name specified."
        return

    if shape_to.observed is None :
        print "Error in plotMappedShapes() : no observed touch points found for {0}. Please check data files.".format(shape_to.name)
        return

    print "Plotting {0} (touch points from {1})".format(shape_to.name, shape_to.pair_mapping)

    out_path = os.path.join(out_path, patient, "figures", "pair_maps")
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    edge_points_to = shape_to.edge_points
    medial_axis_to = shape_to.medial_axis

    shape_from = Shape(shape_to.img_path, shape_to.pair_mapping)

    edge_points_from = shape_from.edge_points
    medial_axis_from = shape_from.medial_axis

    tf = shape_to.fitted_transforms[shape_from.name]
    ma_rot = rotate_points(medial_axis_from, tf, shape_from.dims)
    ep_rot = rotate_points(edge_points_from, tf, shape_from.dims)

    ob_rot = shape_to.observed

    img_bounds = np.array( [ [0,0], [shape_to.dims[0], 0] , [shape_to.dims[0],shape_to.dims[1]], [0,shape_to.dims[1]] ])
    circ_cent, circ_rad = cv2.minEnclosingCircle(img_bounds)
    circ_bd_x = [np.floor(circ_cent[1] - circ_rad), np.ceil(circ_cent[1] + circ_rad)]
    circ_bd_y = [np.floor(circ_cent[0] - circ_rad), np.ceil(circ_cent[0] + circ_rad)]
    canvas_size = circ_bd_x[1] - circ_bd_x[0]
    canvas_offset = np.array([int(circ_bd_y[0]),int(circ_bd_x[0])]) #(x,y)
    
    canvas = np.zeros((int(canvas_size),int(canvas_size),4), dtype=np.uint8)
    
    circ_cent = np.array([int(circ_cent[0]),int(circ_cent[1])])
    cv2.circle( canvas, tuple(circ_cent-canvas_offset), int(circ_rad), (255,255,255,255), thickness=-1 )
    
    for ep in edge_points_to :
        e = ep - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( canvas, (int(e[0]),int(e[1])), 2, (0,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
    
    # Draw medial axis of shape_to
    for mp in medial_axis_to :
        m = mp - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( canvas, (int(m[0]),int(m[1])), 2, (255,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
        
    # Draw transformed medial axis and edge of shape_from
    for mp in ma_rot :
        m = mp - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( canvas, (int(m[0]),int(m[1])), 2, (0,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
        
    for ep in ep_rot :
        e = ep - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( canvas, (int(e[0]),int(e[1])), 2, (0,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
        
     # Transform observed touchpoints from shape_from to shape_to   
    for op in ob_rot :
        o = op - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( canvas, (int(o[0]),int(o[1])), 2, (0,255,0,255), thickness=-1, lineType=cv2.LINE_AA )
    
    out_fname = "_".join((shape_to.pair_mapping, "to", shape_to.name, "Patient", patient, 'pair_map.png'))

    cv2.imwrite( os.path.join(out_path,out_fname), canvas )


if __name__ == '__main__':

    img_names  = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
                  "blake_01","blake_03","blake_04","blake_06","blake_07",
                  "blake_08","blake_09","blake_10","blake_11","blake_12"]

    patients   = ["DF","MC"]

    in_path  = "D:/ShapeTapper-Analysis/"
    out_path = "D:/ShapeTapper-Analysis/"

    shapes = ShapeIO(in_path, out_path, img_names)

    shapes.run(plotHeatMap, patients)

    shapes.run(plotMedialAxis)