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
import distances as dist
from timeit import default_timer as timer

################## Globals - CHANGE THESE TO RUN ON SPECIFIC SUBJECTS AND SHAPE SETS
analysis_conds = ["bounding_circle","in_shape","touchpoint_hull","patient_fitted"]
img_names = ["blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12",
             "solo5","solo6","solo7","solo9","solo10","solo11","solo12"]
patients = ["DF","MC","MC2"]   
img_path = "./Shapes/"         # path containing shape images
out_path_prefix = "D:/ShapeTapper-Analysis/"

# Parameters
n_sets = 100000   # number of uniform data sets to generate; each data set contains a number of uniform points
                  # equal to the number of observed points for that shape



################# Function Definitions #########################################################################
def load_mat(mat_path, mat_name, img_name) :
    try:
        dat = sio.loadmat(os.path.join(mat_path, mat_name))
    except (TypeError, IOError) :
        print "Error loading: {0} -- skipping {1}...\n".format(mat_name, img_name)
        return None
    return dat


def plot_generated(generated, img) :
    for i in range(generated.shape[0]):
        for p in generated[i] :
            cv2.circle( img, tuple([int(p[0]),int(p[1])]), 1, (0,0,255,255), thickness=-1 ) 
    cv2.namedWindow('image', cv2.WINDOW_NORMAL)
    cv2.imshow('image', img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()
    

def get_scaling_factor(observed, edge_points, centroid) :
    # Retrieve set of all touch points out of shape bounds
    observed_oob = []
    for i, tp in enumerate(observed) :
        if cv2.pointPolygonTest(edge_points, tuple(tp), measureDist=False) == -1 :
            observed_oob.append(tp)
    observed_oob = np.array(observed_oob, dtype=np.float32)
    
    if observed_oob.shape[0] == 0 :
        return 1
    else :
        # Percentage of total touch points which are oob
        percent_oob = np.true_divide(observed_oob.shape[0], observed.shape[0])

        # Percentile reference of OOB distances
        dist_idx    = int(np.ceil(percent_oob * observed_oob.shape[0]))
        if observed_oob.shape[0] == 1 and dist_idx == 1 : dist_idx = 0

        # Find distance corresponding to percentile reference of all OOB distances from the shape's edge
        oob_dists   = dist.points2points(observed_oob, edge_points)
        dilate_dist = np.sort(oob_dists)[dist_idx]

        # Find edge point closest to (dist_idx)-th most OOB point
        edge_pt_ref = np.where(dist.points2points(edge_points,observed_oob) == dilate_dist)

        # May be necessary to increase the interval of equality due to float rounding error
        almost_equal = 0.001 # increase gradually to avoid finding multiple edge points
        while edge_pt_ref[0].shape[0] == 0 :
            edge_pt_ref = np.where(np.abs(dist.points2points(edge_points,observed_oob)-dilate_dist)<=almost_equal)
            almost_equal += 0.001

        # Scaling factor is the amount which that edge point must be expanded to fit (dist_idx)-th percentile point
        edge_pt_ref_dist = dist.points2points(np.array([edge_points[edge_pt_ref]]),centroid)
        return (edge_pt_ref_dist[0] + dilate_dist) / edge_pt_ref_dist[0]


def dilate_edges(edge_points, centroid, scale) :
   edges_centered = edge_points - centroid
   edges_scaled = (edges_centered * scale).astype(int) # floor
   return edges_scaled + centroid


# generates uniformly distributed data points inside a region defined by edge_points
def gen_uniform_points_bounds(n_sets, n_pts, edge_points, mins, maxes) :
    uniform_points = np.zeros((n_sets,n_pts,2)) # output array

    for i in range(n_sets) :
        dist.gen_uniform_bounds(edge_points, mins[1], maxes[1], mins[0], maxes[0], n_pts, uniform_points[i])
    return uniform_points


# generates uniformly distributed data points inside a circle
def gen_uniform_points_circle(n_sets, n_pts, center, radius):
    uniform_points = np.empty((n_sets,n_pts,2)) # output array
    
    for i in range(n_sets) :
        dist.gen_uniform_circle(center[1], center[0], radius, n_pts, uniform_points[i])

    return uniform_points


def gen_points(img_name, img_path, img_mat, mat_path, obs_mat, obs_path, out_path) :
    ################### Load data ###################

    # Read in the image 
    img_file = img_name + ".png"
    img = cv2.imread(os.path.join(img_path, img_file),cv2.IMREAD_UNCHANGED)
    img[(img[:,:,3]==0),0:3] = 0 # Convert alpha transparency to black

    # Load shape analysis data (medial axis, edge, centroid)
    shape_analysis = load_mat(mat_path, img_mat, img_name)
    if (shape_analysis is None) : return
    edge_points = np.ascontiguousarray(shape_analysis['edge_points'], dtype=np.float32)
    centroid = shape_analysis['centroid'].astype(np.float32)
    
    # Load observed touchpoint data
    observed_mat = load_mat(obs_path, obs_mat, img_name)
    if (observed_mat is None) : return
    observed = np.ascontiguousarray(observed_mat['img_dataset'], dtype=np.float32)
    if observed.shape[0] == 0 : return
    observed[:,1] = img.shape[0]-observed[:,1] # opencv coordinates use inverted y-axis
    

    ################### Random point generation ################### 

    generated_data_sets = []

    # Generate uniform data in circle around bounding rectangle (same as touchpoint roi)
    if cond == "bounding_circle" :
        start = timer()
        n_pts = observed.shape[0]
        img_bounds = np.array( [ [0,0], [img.shape[0], 0] , [img.shape[0],img.shape[1]], [0,img.shape[1]] ])
        bd_circ_cent, bd_circ_rad = cv2.minEnclosingCircle(img_bounds)
        generated_data_sets = gen_uniform_points_circle(n_sets, n_pts, bd_circ_cent, bd_circ_rad)

    # Generate uniform data within shape   
    if cond == "in_shape" :
        start = timer()
        # Keep only touch points inside of shape
        observed_inshape = []
        for i, tp in enumerate(observed) :
            if cv2.pointPolygonTest(edge_points, tuple(tp), measureDist=False) == 1 :
                observed_inshape.append(tp)
        observed = np.array(observed_inshape, dtype=np.float32) 
        n_pts = observed.shape[0]
        generated_data_sets = gen_uniform_points_bounds(n_sets, n_pts, edge_points, (0,0), img.shape)

    # Generate uniform data within convex hull of touchpoints    
    if cond == "touchpoint_hull" :
        start = timer()
        n_pts = observed.shape[0]
        observed_hull = np.squeeze(cv2.convexHull(observed))
        hull_min = tuple( [np.min(observed_hull[:,1]), np.min(observed_hull[:,0])] )
        hull_max = tuple( [np.max(observed_hull[:,1]), np.max(observed_hull[:,0])] )
        generated_data_sets = gen_uniform_points_bounds(n_sets, n_pts, observed_hull, hull_min, hull_max)

    # Generate uniform data within area defined by percentage and distance of touchpoints outside shape
    if cond == "patient_fitted" :
        start = timer()
        n_pts = observed.shape[0]

        # Dialate edge points to include points within a certain limit outside the shape
        scale = get_scaling_factor(observed, edge_points, centroid)
        edge_points_d = dilate_edges(edge_points, centroid, scale)

        mins_d  = tuple( [centroid[0,1] - img.shape[0]*scale/2, centroid[0,0] - img.shape[1]*scale/2] )
        maxes_d = tuple( [centroid[0,1] + img.shape[0]*scale/2, centroid[0,0] + img.shape[1]*scale/2] )
        generated_data_sets = gen_uniform_points_bounds(n_sets, n_pts, edge_points_d, mins_d, maxes_d)

    # Sanity checks
    print "Generated", img_name, "in", timer()-start, "s", generated_data_sets.shape
    # plot_generated(generated_data_sets, img)

    sio.savemat(out_path+img_name+'_Patient_'+patient+'_generated_uniform_sets.mat', 
                {'gen_datasets':generated_data_sets})



if __name__ == '__main__':              
    for patient in patients :
        print '\n', "Patient:", patient

        mat_path = img_path+"shape_analysis/"
        obs_path = "./"+patient+"/aggregated_observations/"

        for cond in analysis_conds :
            print '\n', "Condition:", cond

            out_path = os.path.join(out_path_prefix,patient+"/generated_uniform_data/"+cond+"/")
            try:
                os.makedirs(out_path)
            except OSError as e:
                if e.errno != errno.EEXIST:
                    raise

            for img_name in img_names :
                img_file = img_name + ".png"
                img_mat  = img_name + "_shape_analysis.mat"
                obs_mat  = img_name + "_Patient_"+patient+"_aggregated_observations.mat"

                gen_points(img_name, img_path, img_mat, mat_path, obs_mat, obs_path, out_path)