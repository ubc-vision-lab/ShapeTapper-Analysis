#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 20 14:50:39 2018

Generates sets of uniform data points corresponding to three shape masks:

Bounding circle - All points within the circle defined by the diagonal of the
                  shape's bounding rectangle (image canvas)
In shape        - All points within the boundary of the shape itself
Touchpoint hull - All points within a space defined by the perimeter of the touchpoints

These sets are specific to each patient's data set, so this script must be run before
running mat_analysis.py

Generated sets will be saved in the directory "generated_uniform_sets" with a subdirectory
for each shape mask

@author: Jamie Dunkle
"""

import cv2
import os, errno
import numpy as np
import scipy.io as sio
from _get_dists import ffi, lib
from timeit import default_timer as timer

################## Globals - CHANGE THESE TO RUN ON SPECIFIC SUBJECTS AND SHAPE SETS
analysis_conds = ["bounding_circle","in_shape","touchpoint_hull"]
img_names = ["blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12",
            "solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12"]
patient = "DF"   
img_path = "./Shapes/"         # path containing shape images

# Parameters
n_sets = 100000   # number of uniform data sets to generate; each data set contains a number of uniform points
                  # equal to the number of observed points for that shape
#scale = 1.2      # expansion factor to dilate shape (CURRENTLY UNUSED)



################# Function Definitions #########################################################################
def plot_generated(generated, img) :
    for i in range(generated.shape[0]):
        for p in generated[i] :
            cv2.circle( img, tuple([int(p[0]),int(p[1])]), 1, (0,0,255,255), thickness=-1 ) 
    cv2.namedWindow('image', cv2.WINDOW_NORMAL)
    cv2.imshow('image', img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    
    
#def dilateEdges(edge_points, centroid, scale) :
#    edges_centered = edge_points - centroid
#    edges_scaled = (edges_centered * scale).astype(int) # floor
#    return edges_scaled + centroid


# generates uniformly distributed data points inside a region defined by edge_points
def gen_uniform_points_bounds(n_sets, n_pts, edge_points, img_dims_min, img_dims_max) :
    uniform_points = np.zeros((n_sets,n_pts,2)) # output array

    # C subroutine uses random seeds, and is called too frequently for system time seeds
    seeds = np.random.random_integers(np.iinfo(np.int32).min, np.iinfo(np.int32).max, n_sets)
    e_pts = ffi.cast("float *", edge_points.ctypes.data) # random set
    rand_decl = "float[{0}]".format(n_pts*2) # string to allocate float array of constant size
    
    for i in range(n_sets) :
        rands = ffi.new(rand_decl) # random set
        lib.gen_uniform_bounds(rands, e_pts, edge_points.shape[0], img_dims_min[1], img_dims_max[1], img_dims_min[0], img_dims_max[0], n_pts, seeds[i])
        uniform_points[i] = np.array(ffi.unpack(rands, n_pts*2), dtype=np.float32).reshape((n_pts,2))

    return uniform_points


# generates uniformly distributed data points inside a circle
def gen_uniform_points_circle(n_sets, n_pts, circ_cent, circ_rad):
    uniform_points = np.zeros((n_sets,n_pts,2)) # output array

    seeds = np.random.random_integers(np.iinfo(np.int32).min, np.iinfo(np.int32).max, n_sets)
    rand_decl = "float[{0}]".format(n_pts*2) # string to allocate float array of constant size

    for i in range(n_sets) :
        rands = ffi.new(rand_decl) # random set
        lib.gen_uniform_circle(rands, circ_cent[1], circ_cent[0], circ_rad, n_pts, seeds[i])
        uniform_points[i] = np.array(ffi.unpack(rands, n_pts*2), dtype=np.float32).reshape((n_pts,2))

    return uniform_points



if __name__ == '__main__':              

    mat_path = img_path+"shape_analysis/"
    dat_path = "./"+patient+"/aggregated_observations/"

    for cond in analysis_conds :
        out_path = "./"+patient+"/generated_uniform_data/"+cond+"/"
        try:
            os.makedirs(out_path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        print "Condition:", cond
        
        for img_name in img_names :
            generated_data_sets = []
            
            img_file = img_name + ".png"
            img_mat = img_name + "_shape_analysis.mat"
            data_mat = img_name + "_Patient_"+patient+"_aggregated_observations.mat"
    
            try:
                shape_analysis = sio.loadmat(os.path.join(mat_path, img_mat))
            except (TypeError, IOError) :
                print "Error loading: {0} -- skipping {1}...\n".format(img_mat, img_name)
                continue
            edge_points = np.ascontiguousarray(shape_analysis['edge_points'], dtype=np.float32)
            # centroid = shape_analysis['centroid'].astype(np.int32)
            
            try:
                observed_mat = sio.loadmat(os.path.join(dat_path, data_mat))
            except (TypeError, IOError) :
                print "Error loading: {0} -- skipping {1}...\n".format(data_mat, img_name)
                continue
            observed = np.ascontiguousarray(observed_mat['img_dataset'], dtype=np.float32)
            if observed.shape[0] == 0 : continue

            img = cv2.imread(os.path.join(img_path, img_file),cv2.IMREAD_UNCHANGED)
            img[(img[:,:,3]==0),0:3] = 0 # Convert alpha transparency to black
          
    #        # Dialate edge points to include points within a certain limit outside the shape
    #        edge_points_d = dilateEdges(edge_points, centroid, scale)
    #        img_dims_d = (int(img.shape[0] * scale/2), int(img.shape[1] * scale/2), img.shape[2])

            observed[:,1] = img.shape[0]-observed[:,1] # opencv coordinates use inverted y-axis
            
            if cond == "bounding_circle" :
                start = timer()
                # Generate uniform data in circle around bounding rectangle (same as touchpoint roi)
                n_pts = observed.shape[0]
                img_bounds = np.array( [ [0,0], [img.shape[0], 0] , [img.shape[0],img.shape[1]], [0,img.shape[1]] ])
                bd_circ_cent, bd_circ_rad = cv2.minEnclosingCircle(img_bounds)
                generated_data_sets = gen_uniform_points_circle(n_sets, n_pts, bd_circ_cent, bd_circ_rad)
                
            if cond == "in_shape" :
                start = timer()
                # Generate uniform data within shape
                observed_inshape = []
                for i in range(observed.shape[0]) :
                    if cv2.pointPolygonTest(edge_points, tuple(observed[i]), measureDist=False) == 1 :
                        observed_inshape.append(observed[i])
                observed = np.array(observed_inshape, dtype=np.float32) 
                n_pts = observed.shape[0]
                generated_data_sets = gen_uniform_points_bounds(n_sets, n_pts, edge_points, (0,0), img.shape)
                
            if cond == "touchpoint_hull" :
                start = timer()
                # Generate uniform data within convex hull of touchpoints
                n_pts = observed.shape[0]
                observed_hull = np.squeeze(cv2.convexHull(observed))
                hull_min = tuple( [np.min(observed_hull[:,1]), np.min(observed_hull[:,0])] )
                hull_max = tuple( [np.max(observed_hull[:,1]), np.max(observed_hull[:,0])] )
                generated_data_sets = gen_uniform_points_bounds(n_sets, n_pts, observed_hull, hull_min, hull_max)

            sio.savemat(out_path+img_name+'_Patient_'+patient+'_generated_uniform_sets.mat', {'gen_datasets':generated_data_sets})

            # Sanity checks
            print "Generated", img_name, "in", timer()-start, "s", generated_data_sets.shape
            # plot_generated(generated_data_sets, img)