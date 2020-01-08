#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Generates plots of touch points, medial axis and analysis data for each shape

@author: Jamie Dunkle
"""

import cv2
import os, errno
import numpy as np
from ShapeObj import Shape
from ShapeIO import ShapeIO
from map_shape_pairs import rotate_points


################# Function Definitions #########################################################################

# Plot color-graded "heat map" representing the relative density of touch points in a given shape
def plotHeatMap(shape, out_path, patient, cond=None, task=None) :

    print "Plotting {0}".format(shape.name)

    if patient is None :
        print "Error in plotHeatMap() : no patient name specified."
        return

    if shape.observed is None :
        print "Error in plotHeatMap() : no observed touch points found for {0}. Please check data files.".format(shape.name)
        return

    # out_path = os.path.join(out_path, patient, "figures", "heat_maps")
    out_path = os.path.join(out_path, "figures", "heat_maps", patient)
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # Fetch observed touch points and Reference Object points from shape 
    observed = shape.observed
    medial_axis = shape.medial_axis.astype(np.uint32)
    edge_points = shape.edge_points.astype(np.uint32)
    centroid = shape.centroid.astype(np.uint32)

    # Create a blank greyscale canvas to contain the bonding circle of the image canvas 
    img_bounds = np.array( [ [0,0], [shape.dims[0], 0] , [shape.dims[0],shape.dims[1]], [0,shape.dims[1]] ])
    circ_cent, circ_rad = cv2.minEnclosingCircle(img_bounds)
    circ_rad = circ_rad * 1.2
    circ_bd_x = [np.floor(circ_cent[1] - circ_rad), np.ceil(circ_cent[1] + circ_rad)]
    circ_bd_y = [np.floor(circ_cent[0] - circ_rad), np.ceil(circ_cent[0] + circ_rad)]
    canvas_size = circ_bd_x[1] - circ_bd_x[0]
    canvas_offset = np.array([int(circ_bd_y[0]),int(circ_bd_x[0])]) #(x,y)
    canvas = np.zeros((int(canvas_size),int(canvas_size)), dtype=np.uint8)
    
    # Draw black circle to represent bounding circle
    circ_cent = np.array([int(circ_cent[0]),int(circ_cent[1])])
    cv2.circle( canvas, tuple(circ_cent-canvas_offset), int(circ_rad), 0, thickness=-1 )
    
    # Plot observed touch points with alpha-channel fade
    for op in observed :
        o = np.array([int(op[0]),int(op[1])]) - np.flip(canvas_offset,axis=0) #(y,x)
        overlay = canvas.copy()
        rad = (int) (np.min(shape.dims[0:2]) / 12)#(observed.shape[0]/22)) # point radius
        for i in range(1, rad+1) :
            alpha = np.true_divide(1,4*rad) # alpha decreases with radius
            cv2.circle( overlay, tuple(o), i, 255, thickness=-1, lineType=cv2.LINE_AA )
            cv2.addWeighted(overlay, alpha, canvas, 1 - alpha, 0, canvas) # add image to canvas with alpha-channel
    
    # Normalize touchpoints to full range of grey values (to ensure each shape has full color spectrum)
    canvas_normed = canvas.copy()
    cv2.normalize(canvas, canvas_normed, alpha=0, beta=255, norm_type=cv2.NORM_MINMAX)

    # Apply color map ("JET") to greyscale image to produce heat_map
    heat_map = np.zeros((int(canvas_size),int(canvas_size),4), dtype=np.uint8)
    heat_map = cv2.applyColorMap(canvas_normed, cv2.COLORMAP_JET)

    # Create add an alpha channel to heat_map
    b_channel, g_channel, r_channel = cv2.split(heat_map)
    alpha_channel = np.ones(b_channel.shape, dtype=b_channel.dtype) * 255 
    heat_map = cv2.merge((b_channel, g_channel, r_channel, alpha_channel))
    
    # Get blue pixels (bottom of JET color map) & convert to transparent
    heat_map[heat_map[:,:,0]==128] = [0,0,0,0] 

    # Use bitwise masks to copy heat_map onto the bounding circle canvas, with transparent background 
    img_mask = np.zeros((int(canvas_size),int(canvas_size),4), dtype=np.uint8)
    cv2.fillConvexPoly(img_mask, edge_points - np.flip(canvas_offset,axis=0) , (255,255,255,255))
    img_mask = cv2.bitwise_or(img_mask,heat_map)
    heat_map[heat_map[:,:,3]==0] = [255,255,255,255]
    heat_map = cv2.bitwise_and(heat_map,img_mask)

    # Set line width (for edges) according to image dimensions
    line_width = np.max(shape.dims[0:2])/95 

    # Plot medial axis, centroid and edges
    for mp in medial_axis :
        m = mp - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( heat_map, tuple(m), line_width, (255,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
    
    c = centroid[0] - np.flip(canvas_offset,axis=0) #(y,x)
    cv2.circle( heat_map, tuple(c), line_width*3, (0,255,0,255), thickness=-1, lineType=cv2.LINE_AA )  

    # heat_map[heat_map[:,:,3]==0] = [0,0,0,255] # convert shape interior to transparent

    for ep in edge_points :
        e = ep - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( heat_map, tuple(e), line_width, (0,0,0,255), thickness=-1, lineType=cv2.LINE_AA )

    # Generate output file name
    if task is None:
        if shape.pair_mapping is not None :
            out_fname = "_".join((shape.pair_mapping, "to", shape.name, "Patient", patient, 'heat_map.png'))
        else :
            out_fname = "_".join((shape.name, "Patient", patient, 'heat_map.png'))
    else :
        if shape.pair_mapping is not None :
            out_fname = "_".join((shape.pair_mapping, "to", shape.name, "Patient", patient, task, '_heat_map.png'))
        else :
            out_fname = "_".join((shape.name, "Patient", patient, task, 'heat_map.png'))

    # Save heat_map as PNG
    cv2.imwrite( os.path.join(out_path,out_fname), heat_map )


# Plots the shape with its medial axis and centroid
def plotRefObjects(shape, out_path, patient=None, cond=None, task=None) :

    print "Plotting {0}".format(shape.name)

    out_path = os.path.join(out_path, "Shapes", "figures", "medial_axis")
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # Fetch observed touch points and Reference Object points from shape 
    img = shape.img.copy()
    medial_axis = shape.medial_axis.astype(np.uint32)
    edge_points = shape.edge_points.astype(np.uint32)
    centroid = shape.centroid.astype(np.uint32)
    
    # Set line width according to image canvas dimensions
    line_width = np.max(img.shape[0:2])/190 

    # Plot RefObj points
    for mp in medial_axis :
       cv2.circle( img, tuple(mp), line_width, (255,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
   
    cv2.circle( img, tuple(centroid[0]), line_width*4, (0,255,0,255), thickness=-1, lineType=cv2.LINE_AA )  
    
    for ep in edge_points :
        cv2.circle( img, tuple(ep), line_width, (0,0,0,255), thickness=-1, lineType=cv2.LINE_AA )

    # heat_map[heat_map[:,:,3]==0] = [0,0,0,255] # convert shape interior to transparent

    # Generate output file name and save image as PNG
    out_fname = "_".join((shape.name, '_ma.png'))
    cv2.imwrite( os.path.join(out_path, out_fname), img )


# Plots the mapped shape pairs produced by map_shape_pairs
def plotMappedShapes(shape_to, out_path, patient, cond=None, task=None) :

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

    # Fetch observed touch points and Reference Object points from shape 
    edge_points_to = shape_to.edge_points
    medial_axis_to = shape_to.medial_axis

    # Load shape_from and fetch Reference Object points
    shape_from = Shape(shape_to.img_path, shape_to.pair_mapping)
    edge_points_from = shape_from.edge_points
    medial_axis_from = shape_from.medial_axis

    # Fetch transform params to shape_from and rotate points to map to shape_to
    tf = shape_to.fitted_transforms[shape_from.name]
    ma_rot = rotate_points(medial_axis_from, tf, shape_from.dims)
    ep_rot = rotate_points(edge_points_from, tf, shape_from.dims)

    # Transformed observed touch points will already be loaded into the shape object by ShapeIO
    ob_rot = shape_to.observed

    # Create a blank canvas to contain the bonding circle of the image canvas 
    img_bounds = np.array( [ [0,0], [shape_to.dims[0], 0] , [shape_to.dims[0],shape_to.dims[1]], [0,shape_to.dims[1]] ])
    circ_cent, circ_rad = cv2.minEnclosingCircle(img_bounds)
    circ_bd_x = [np.floor(circ_cent[1] - circ_rad), np.ceil(circ_cent[1] + circ_rad)]
    circ_bd_y = [np.floor(circ_cent[0] - circ_rad), np.ceil(circ_cent[0] + circ_rad)]
    canvas_size = circ_bd_x[1] - circ_bd_x[0]
    canvas_offset = np.array([int(circ_bd_y[0]),int(circ_bd_x[0])]) #(x,y)
    
    canvas = np.zeros((int(canvas_size),int(canvas_size),4), dtype=np.uint8)
    
    # Draw white circle to represent bounding circle
    circ_cent = np.array([int(circ_cent[0]),int(circ_cent[1])])
    cv2.circle( canvas, tuple(circ_cent-canvas_offset), int(circ_rad), (255,255,255,255), thickness=-1 )
    
    # Draw edges of shape_to
    for ep in edge_points_to :
        e = ep - np.flip(canvas_offset,axis=0) # (y,x), offset to canvas
        cv2.circle( canvas, (int(e[0]),int(e[1])), 2, (0,0,0,255), thickness=-1, lineType=cv2.LINE_AA )
    
    # Draw medial axis of shape_to
    for mp in medial_axis_to :
        m = mp - np.flip(canvas_offset,axis=0) # (y,x), offset to canvas
        cv2.circle( canvas, (int(m[0]),int(m[1])), 2, (255,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
        
    # Draw transformed medial axis and edge of shape_from
    for mp in ma_rot :
        m = mp - np.flip(canvas_offset,axis=0) # (y,x), offset to canvas
        cv2.circle( canvas, (int(m[0]),int(m[1])), 2, (255,0,0,255), thickness=-1, lineType=cv2.LINE_AA )
    for ep in ep_rot :
        e = ep - np.flip(canvas_offset,axis=0) # (y,x), offset to canvas
        cv2.circle( canvas, (int(e[0]),int(e[1])), 2, (0,0,255,255), thickness=-1, lineType=cv2.LINE_AA )
        
     # Transform observed touchpoints from shape_from to shape_to   
    for op in ob_rot :
        o = op - np.flip(canvas_offset,axis=0) # (y,x), offset to canvas
        cv2.circle( canvas, (int(o[0]),int(o[1])), 2, (0,255,0,255), thickness=-1, lineType=cv2.LINE_AA )
    
    # Generate output name and save figure as PNG
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

    shapes.run(plotRefObjects)