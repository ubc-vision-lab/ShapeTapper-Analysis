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
import scipy.io as sio


################## Globals - CHANGE THESE TO RUN ON SPECIFIC SUBJECTS AND SHAPE SETS
img_names = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
                "blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"]
patient = "MC"
img_path = "./Shapes/" 

# Heatmap parameters
rad = 12                    # Touchpoint radius in pixels (default = 12)
alpha = 0.033               # Transparency value (default=0.033)
# h_color = (0,0,255,255)   # Touchpoint color (default=(0,155,255,255), orange)


################# Function Definitions #########################################################################
def makeHeatMap(img_dims, edge_points, ma_points, centroid, observed) :
    
    img_bounds = np.array( [ [0,0], [img_dims[0], 0] , [img_dims[0],img_dims[1]], [0,img_dims[1]] ])
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
        rad = (int) (np.min(img_dims[0:2]) / 12)
        alpha = np.true_divide(1,2*rad)
        for i in range(1, rad+1) :
            cv2.circle( overlay, tuple(o), i, 255, thickness=-1, lineType=cv2.LINE_AA )
            cv2.addWeighted(overlay, alpha, canvas, 1 - alpha, 0, canvas)  

    heat_map = cv2.applyColorMap(canvas, cv2.COLORMAP_JET)

    # heat_map[heat_map[:,:,0]==128] = [255,255,255]

    for ep in edge_points :
        e = ep - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( heat_map, tuple(e), 2, (255,255,255,255), thickness=-1, lineType=cv2.LINE_AA )
    
    for mp in ma_points :
        m = mp - np.flip(canvas_offset,axis=0) #(y,x)
        cv2.circle( heat_map, tuple(m), 2, (255,0,255,255), thickness=-1, lineType=cv2.LINE_AA )

    c = centroid[0] - np.flip(canvas_offset,axis=0) #(y,x)
    cv2.circle( heat_map, tuple(c), 6, (0,255,0,255), thickness=-1, lineType=cv2.LINE_AA )    

    # cv2.namedWindow('image', cv2.WINDOW_NORMAL)
    # cv2.imshow('image', heat_map)
    # cv2.waitKey(0)
    # cv2.destroyAllWindows()
    return heat_map



if __name__ == '__main__':

    mat_path = img_path+"shape_analysis/"               # path containing medial axis mat files
    obs_path = "./"+patient+"/aggregated_observations/" # path containing observed data
    
    out_path = "./"+patient+"/heat_maps/"       # output path
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    for img_name in img_names :
        
        print 'Starting', img_name
        img_file = img_name + ".png"
        img_mat = img_name + "_shape_analysis.mat"
        obs_mat = img_name + "_Patient_"+patient+"_aggregated_observations.mat"
        
        img = cv2.imread(os.path.join(img_path, img_file) ,cv2.IMREAD_UNCHANGED)
        img[(img[:,:,3]==0),0:3] = 0
        
        try:
            shape_analysis = sio.loadmat(os.path.join(mat_path, img_mat))
        except (TypeError, IOError) :
            print "Error loading: {0} -- skipping {1}...".format(img_mat, img_name)
            continue
        ma_points = shape_analysis['ma_points'].astype(np.int32) # (x,y)
        edge_points = shape_analysis['edge_points'].astype(np.int32)  # (x,y)
        centroid = shape_analysis['centroid'].astype(np.int32)  # (x,y)
        
        try:
            observed_mat = sio.loadmat(os.path.join(obs_path, obs_mat))
        except (TypeError, IOError) :
            print "Error loading: {0} -- skipping {1}...".format(obs_mat, img_name)
            continue
        observed = observed_mat['img_dataset'].astype(np.float32)
        observed[:,1] = img.shape[0]-observed[:,1] # opencv coordinates use inverted y-axis
       
        heat_map = makeHeatMap(img.shape, edge_points, ma_points, centroid, observed)
        cv2.imwrite(os.path.join(out_path, img_name + '_heat_map.png'),heat_map)