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


################## Globals - CHANGE THESE TO RUN ON SPECIFIC SUBJECTS AND SHAPE SETS
analysis_conds = ["bounding_circle","in_shape","touchpoint_hull"]
img_names = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
             "blake_01","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"]
patient = "MC"   
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
# using the include/exclude method via pointPolygonTest
def gen_uniform_points(n_sets, n_pts, edge_points, img_dims, offset = 10):
    num_uniform = np.zeros(n_sets, dtype=np.uint8) # count num of included points
    uniform_points = np.zeros((n_sets,n_pts,2))
    edge_points += offset # pointPolygonTest fails if edge points fall on img canvas borders
    while True:
        # generate random data points in the shape's bounding rectangle
        rand_set_y = np.random.uniform (0, img_dims[0]+1, (n_pts*2,n_sets))
        rand_set_x = np.random.uniform (0, img_dims[1]+1, (n_pts*2,n_sets))
        rand_points = np.array([rand_set_x, rand_set_y]).transpose()
        rand_points += offset
        for i in range(rand_points.shape[0]) : # for each point set
            for rpt in rand_points[i] :        # for each point
                if num_uniform[i] >= n_pts :
                        continue # skip any point set which has been filled
                else:
                    if cv2.pointPolygonTest(edge_points, tuple(rpt), measureDist=False) == 1 :
                        uniform_points[i,num_uniform[i],:] = np.array(rpt)
                        num_uniform[i] += 1
        if np.all(num_uniform >= n_pts): # if not all sets have been filled, loop again
            uniform_points -= offset
            return uniform_points
        
        
# generates uniformly distributed data points inside a circular region
# using the include/exclude method via pointPolygonTest
def gen_uniform_points_circ(n_sets, n_pts, circ_cent, circ_rad):
    num_uniform = np.zeros(n_sets, dtype=np.uint8) # count num of included points
    uniform_points = np.zeros((n_sets,n_pts,2))
    circ_bd_x = [np.floor(circ_cent[1] - circ_rad), np.ceil(circ_cent[1] + circ_rad)]
    circ_bd_y = [np.floor(circ_cent[0] - circ_rad), np.ceil(circ_cent[0] + circ_rad)]
    while True:
        # generate random data points in the shape's bounding rectangle
        rand_set_y = np.random.uniform (circ_bd_y[0], circ_bd_y[1], (n_pts*2,n_sets))
        rand_set_x = np.random.uniform (circ_bd_x[0], circ_bd_x[1], (n_pts*2,n_sets))
        rand_points = np.array([rand_set_y, rand_set_x]).transpose()
        for i in range(rand_points.shape[0]) : # for each point set
            for rpt in rand_points[i] :        # for each point
                if num_uniform[i] >= n_pts :
                        continue # skip any point set which has been filled
                else:
                    if np.linalg.norm(rpt - circ_cent) <= circ_rad :
                        uniform_points[i,num_uniform[i],:] = np.array(rpt)
                        num_uniform[i] += 1
        if np.all(num_uniform >= n_pts): # if not all sets have been filled, loop again
            return uniform_points


if __name__ == '__main__':              
    
    mat_path = "./"+patient+"/shape_analysis/" 
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
                print "Error loading: {0} -- skipping {1}...".format(img_mat, img_name)
                continue
            edge_points = shape_analysis['edge_points'].astype(np.int32)
            centroid = shape_analysis['centroid'].astype(np.int32)
            
            try:
                observed_mat = sio.loadmat(os.path.join(dat_path, data_mat))
            except (TypeError, IOError) :
                print "Error loading: {0} -- skipping {1}...".format(data_mat, img_name)
                continue
            observed = observed_mat['img_dataset'].astype(np.float32)
        
            img = cv2.imread(os.path.join(img_path, img_file),cv2.IMREAD_UNCHANGED)
            img[(img[:,:,3]==0),0:3] = 0 # Convert alpha transparency to black
          
    #        # Dialate edge points to include points within a certain limit outside the shape
    #        edge_points_d = dilateEdges(edge_points, centroid, scale)
    #        img_dims_d = (int(img.shape[0] * scale/2), int(img.shape[1] * scale/2), img.shape[2])

            observed[:,1] = img.shape[0]-observed[:,1] # opencv coordinates use inverted y-axis
            
            if cond == "bounding_circle" :
                # Generate uniform data in circle around bounding rectangle (same as touchpoint roi)
                n_pts = observed.shape[0]
                img_bounds = np.array( [ [0,0], [img.shape[0], 0] , [img.shape[0],img.shape[1]], [0,img.shape[1]] ])
                bd_circ_cent, bd_circ_rad = cv2.minEnclosingCircle(img_bounds)
                generated_data_sets = gen_uniform_points_circ(n_sets, n_pts, bd_circ_cent, bd_circ_rad)
                sio.savemat(out_path+img_name+'_Patient_'+patient+'_generated_uniform_sets.mat', {'gen_datasets':generated_data_sets})
                
            if cond == "in_shape" :
                # Generate uniform data within shape
                observed_inshape = []
                for i in range(observed.shape[0]) :
                    if cv2.pointPolygonTest(edge_points, tuple(observed[i]), measureDist=False) == 1 :
                        observed_inshape.append(observed[i])
                observed = np.array(observed_inshape) 
                n_pts = observed.shape[0]
                generated_data_sets = gen_uniform_points(n_sets, n_pts, edge_points, img.shape)
                sio.savemat(out_path+img_name+'_Patient_'+patient+'_generated_uniform_sets.mat', {'gen_datasets':generated_data_sets})
                
            if cond == "touchpoint_hull" :
                # Generate uniform data within convex hull of touchpoints
                n_pts = observed.shape[0]
                observed_hull = np.squeeze(cv2.convexHull(observed))
                generated_data_sets = gen_uniform_points(n_sets, n_pts, observed_hull, img.shape)
                sio.savemat(out_path+img_name+'_Patient_'+patient+'_generated_uniform_sets.mat', {'gen_datasets':generated_data_sets})
    
            # Sanity checks
            print "Generated", img_name, generated_data_sets.shape
    #        plot_generated(generated_data, img)