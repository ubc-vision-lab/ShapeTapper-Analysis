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
running spatial_analysis.py

Generated sets will be saved in the directory "uniform_points with a subdirectory
for each shape mask

@author: Jamie Dunkle
"""

import cv2
import os, errno
import numpy as np
import scipy.io as sio
import random_point_lib as rp
import distance_lib as dist
from ShapeIO import ShapeIO
from timeit import default_timer as timer


################# Function Definitions #########################################################################
def plot_generated(generated, edge_points, img) :
    for i in range(generated.shape[0]):
        for p in generated[i] :
            cv2.circle( img, tuple([int(p[0]),int(p[1])]), 1, (0,0,255,255), thickness=-1, lineType=cv2.LINE_AA ) 

    for e in edge_points :
        cv2.circle( img, (int(e[0]),int(e[1])), 1, (255,0,0,255), thickness=-1, lineType=cv2.LINE_AA )

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

        # May be necessary to use an interval of equality to find closest edge point, due to float rounding error
        almost_equal = 0.0001 # increase gradually to avoid finding multiple edge points
        while edge_pt_ref[0].shape[0] == 0 :
            edge_pt_ref = np.where(np.abs(dist.points2points(edge_points,observed_oob)-dilate_dist)<=almost_equal)
            almost_equal += 0.0001

        # Scaling factor is the amount which that edge point must be expanded to fit (dist_idx)-th percentile point
        edge_pt_ref_dist = dist.points2points(np.array([edge_points[edge_pt_ref]]),centroid)
        return (edge_pt_ref_dist[0] + dilate_dist) / edge_pt_ref_dist[0]

# Expands edge points
def dilate_edges(edge_points, center_point, scale) :
   edges_centered = edge_points - center_point
   edges_scaled = edges_centered * scale
   return edges_scaled + center_point


# generates uniformly distributed data points inside a region defined by edge_points
def gen_uniform_points_bounds(n_sets, n_pts, edge_points, mins, maxes) :
    uniform_points = np.zeros((n_sets,n_pts,2)) # output array
    for i in range(n_sets) :
        rp.gen_uniform_bounds(edge_points, mins[1], maxes[1], mins[0], maxes[0], n_pts, uniform_points[i])
    return uniform_points


# generates uniformly distributed data points inside a circle
def gen_uniform_points_circle(n_sets, n_pts, center, radius):
    uniform_points = np.empty((n_sets,n_pts,2)) # output array
    for i in range(n_sets) :
        rp.gen_uniform_circle(center[1], center[0], radius, n_pts, uniform_points[i])
    return uniform_points

# generates uniformly distributed data points inside a circle
def gen_uniform_points_cent_normal(n_sets, n_pts, center, bd_center, bd_radius, std) :
    uniform_points = np.empty((n_sets,n_pts,2)) # output array
    for i in range(n_sets) :
        rp.gen_uniform_cent_normal(center[0], center[1], bd_center[1], bd_center[0], bd_radius, n_pts, std[0], std[1], uniform_points[i])
    return uniform_points


################### Random point generation ####################################################################
def generateUniformData(shape, out_path, patient, cond, n_sets = 100000) :

    if patient is None :
        print "Error in generateUniformData() : no patient name specified."
        return

    if cond is None :
        print "Error in generateUniformData() : no condition(s) specified."
        return

    if shape.observed is None :
        print "Error in generateUniformData() : no observed touch points found for {0}. Please check data files.".format(shape.name)
        return

    out_path = os.path.join(out_path, patient, "uniform_points", cond)
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    generated_data_sets = []

    observed = shape.observed
    edge_points = shape.edge_points
    centroid = shape.centroid

    # Generate uniform data in circle around bounding rectangle (same as touchpoint roi)
    if cond == "bounding_circle" :
        start = timer()
        n_pts = observed.shape[0]
        img_bounds = np.array( [ [0,0], [shape.dims[0], 0] , [shape.dims[0],shape.dims[1]], [0,shape.dims[1]] ])
        bd_circ_cent, bd_circ_rad = cv2.minEnclosingCircle(img_bounds)
        generated_data_sets = gen_uniform_points_circle(n_sets, n_pts, bd_circ_cent, bd_circ_rad)

     # Generate data normal distribution from the centroid (within the bounding circle as above)
    if cond == "normal_distribution" :
        start = timer()
        n_pts = observed.shape[0]
        img_bounds = np.array( [ [0,0], [shape.dims[0], 0] , [shape.dims[0],shape.dims[1]], [0,shape.dims[1]] ])
        bd_circ_cent, bd_circ_rad = cv2.minEnclosingCircle(img_bounds)
        observed = observed.astype(np.float64)
        std = np.std(observed, axis=0)
        cent = centroid[0].astype(np.float64)
        generated_data_sets = gen_uniform_points_cent_normal(n_sets, n_pts, cent, bd_circ_cent, bd_circ_rad, std).astype(np.float32)
 
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
        generated_data_sets = gen_uniform_points_bounds(n_sets, n_pts, edge_points, (0,0), shape.dims)

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
        scale  = get_scaling_factor(observed, edge_points, centroid)
        center = np.array([shape.dims[0]/2, shape.dims[1]/2], dtype=np.float32)
        edge_points_d = dilate_edges(edge_points, center, scale)

        mins_d  = tuple( [center[0] - shape.dims[0]*scale/2, center[1] - shape.dims[1]*scale/2] )
        maxes_d = tuple( [center[0] + shape.dims[0]*scale/2, center[1] + shape.dims[1]*scale/2] )
        generated_data_sets = gen_uniform_points_bounds(n_sets, n_pts, edge_points_d, mins_d, maxes_d)

    # Sanity checks
    print "Generated {0} in {1}s {2}".format(shape.name, timer()-start, generated_data_sets.shape)
    # plot_generated(generated_data_sets, edge_points, shape.img.copy())

    if shape.pair_mapping is not None :
        out_fname = "_".join((shape.pair_mapping, "to", shape.name, "Patient", patient, "uniform_points", cond)) + '.mat'
    else :
        out_fname = "_".join((shape.name, "Patient", patient, "uniform_points", cond)) + '.mat'
    
    sio.savemat( os.path.join(out_path, out_fname), {'unif_datasets':generated_data_sets} )



if __name__ == '__main__':

    img_names  = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
                  "blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"]

    conditions = ["bounding_circle","in_shape","touchpoint_hull","patient_fitted"]
    patients   = ["DF","MC"]

    in_path  = "D:/ShapeTapper-Analysis/"
    out_path = "D:/ShapeTapper-Analysis/"

    shapes = ShapeIO(in_path, out_path, img_names)

    shapes.run(generateUniformData, patients, conditions)