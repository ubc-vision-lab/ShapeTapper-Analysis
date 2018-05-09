#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 15:42:30 2018

Performs analysis of distance and point clustering (Sadahiro & Takami, 2001) for each shape.

Distance data will be saved in the directory "distance_analysis" with a subdirectory for each shape mask. 
Each file contains data to create a probability distribution defined by the uniform uniform sets.

The MATLAB script StatAnalysis.m will use these to calculate the significance of the observed data. 

@author: Jamie Dunkle
"""

import cv2
import os, errno
import numpy as np
import scipy.io as sio
from numba import guvectorize
import distance_lib as dist
from ShapeIO import ShapeIO
from timeit import default_timer as timer

################# Function Definitions #########################################################################
def getVarMean(data) :
    a = data.ndim - 1
    var = np.sum(np.power(data,2),axis=a)/(data.shape[a]-1)
    rmse = np.sqrt(var)
    amd = np.mean(data, axis=a)
    return [var, rmse, amd]


def get_dists(points, medial_axis, edge_points, centroid) :
    if points is None : return [[],[],[]]

    if len(points.shape) < 3 :
        ma_dists   = dist.points2points(points,medial_axis)
        edge_dists = dist.points2points(points,edge_points)
        cent_dists = dist.points2points(points,centroid)
    else :
        ma_dists   = np.empty(points.shape[0:2], dtype=np.float32)
        edge_dists = np.empty_like(ma_dists)
        cent_dists = np.empty_like(ma_dists)
        for i, point_set in enumerate(points) :
            ma_dists[i]   = dist.points2points(point_set,medial_axis)
            edge_dists[i] = dist.points2points(point_set,edge_points)
            cent_dists[i] = dist.points2points(point_set,centroid)

    medaxis_data = getVarMean(ma_dists)
    edge_data    = getVarMean(edge_dists)
    cent_data    = getVarMean(cent_dists)

    return [medaxis_data, edge_data, cent_data]


# Cumulative Density Function for Spatial Analysis
# Takes sorted distances for a set of points (ft) & a vector of distances (rs)
# Counts number of points within each distance (rs[r]), as a function cdf
@guvectorize(['void(float32[:], float32[:], float32[:])'],
              "(m),(n)->(m)", nopython=True, cache=True)
def gu_cdf(rs,ft,cdf):
    idx = 0     # counts # of points in given region
    cdf[0] = 0  # counts current region (1:len(rs))
    npts = ft.shape[0]
    for r in range(rs.shape[0]-1) :
        while ft[idx] < rs[r] :
            cdf[r] += 1
            idx += 1    
            if idx == npts :
                cdf[r+1:] = npts # fill remaining CDF with n_pts ==> CDF = 1
                return
        cdf[r+1] = cdf[r] # copy current count to next region

# cumulative distribution function for "reference object" defined by ro_pts
def cdf(points_in, ro_pts, regions):
    ft_dists = np.sort( dist.points2points(points_in,ro_pts) )
    cdf = np.empty(regions.shape[0], dtype=np.float32) # output vector
    gu_cdf(regions, ft_dists, cdf)
    return np.true_divide(cdf, points_in.shape[0]) # normalize [0, n_pts] to [0,1]

def get_cdf(ro_pts, regions, interior_points, observed, uniform=None) :
    # calculate d-values for main shape
    ds_exp = cdf(interior_points, ro_pts, regions)

    # calculate d-values for observed data
    ds_obs = cdf(observed, ro_pts, regions)

    # calculate d-values for 100k simulated data sets
    ds_unif = []
    if uniform is not None :
        ds_unif = np.empty((uniform.shape[0],regions.shape[0]), dtype=np.float32)
        for i, point_set in enumerate(uniform) :
            ds_unif[i] = cdf(point_set, ro_pts, regions)

    return [ds_exp, ds_obs, ds_unif]



def spatialAnalysis(shape, out_path, patient, cond, n_cdf = 1000):    

    if patient is None :
        print "Error in spatialAnalysis() : no patient name specified."
        return

    if cond is None :
        print "Error in spatialAnalysis() : no condition(s) specified."
        return

    if shape.observed is None :
        print "Error in spatialAnalysis() : no observed touch points found for {0}. Please check data files.".format(shape.name)
        return

    if shape.uniform is None :
        print "Warning in spatialAnalysis() : no uniform points found for {0}. Please check data files.".format(shape.name)

    out_path = os.path.join(out_path, patient, "spatial_analysis", cond)
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    observed = shape.observed
    uniform  = shape.uniform

    medial_axis = shape.medial_axis
    edge_points = shape.edge_points
    centroid    = shape.centroid

    # Trim observed touchpoints outside of shape for the "in shape" condition
    if cond == "in_shape" :
        observed_inshape = [] 
        for i in range(observed.shape[0]) :
            if cv2.pointPolygonTest(edge_points, tuple(observed[i]), measureDist=False) == 1 :
                observed_inshape.append(observed[i])
        observed = np.array(observed_inshape)
        
    ################### Spatial Analysis ################### 

    # Calculate variance & mean info for observed and uniform sets
    start = timer()
    observed_varmean = get_dists(observed, medial_axis, edge_points, centroid)
    uniform_varmean  = get_dists(uniform, medial_axis, edge_points, centroid)
    print "{0} analysis : var-mean data generated in {1}s".format(shape.name, timer() - start)
    
    # Find uniformly distributed points in shape (pixel locations)
    img_bin = cv2.cvtColor(shape.img, cv2.COLOR_BGR2GRAY)
    img_bin[(img_bin!=0)] = 1 # convert to logical array (white = 1)
    white_idx = np.where(img_bin==1)    
    shape_points = np.flip(np.transpose(np.array(white_idx)), axis=1)
    shape_points = np.ascontiguousarray(shape_points, dtype=np.float32)

    # Calculate vector of regions for spatial analysis (from 0 to bounds of enclosing circle)
    # Number of steps is specified in NUM_RS variable at top of script
    img_bounds = np.array( [ [0,0], [shape.dims[0], 0] , [shape.dims[0],shape.dims[1]], [0,shape.dims[1]] ])
    _, bd_circ_radius = cv2.minEnclosingCircle(img_bounds)
    max_r_bd = int(np.ceil(bd_circ_radius))
    regions_bd = np.linspace(0.0, max_r_bd, num=n_cdf, endpoint=True)
    regions_bd = np.ascontiguousarray(regions_bd, dtype=np.float32)
    
    # Use a slightly enlarged regions vector for centroid to catch all points
    max_r_ct = int(np.ceil(bd_circ_radius*np.sqrt(2))) 
    regions_ct = np.linspace(0.0, max_r_ct, num=n_cdf, endpoint=True)
    regions_ct = np.ascontiguousarray(regions_ct, dtype=np.float32)

    # Calculate spatial CDF data for each Reference Object (Medax, Edge, Centroid)
    start = timer()
    cdf_ma   = get_cdf(medial_axis, regions_bd, shape_points, observed, uniform)
    cdf_edge = get_cdf(edge_points, regions_bd, shape_points, observed, uniform)
    cdf_cent = get_cdf(centroid, regions_ct, shape_points, observed, uniform)
    print "{0} analysis : cdf data generated in {1}s".format(shape.name, timer() - start)

    # Output data for statistical analysis
    if shape.pair_mapping is not None :
        out_fname = "_".join((shape.pair_mapping, "to", shape.name, "Patient", patient, "spatial_analysis", cond)) + '.mat'
    else :
        out_fname = "_".join((shape.name, "Patient", patient, "spatial_analysis", cond)) + '.mat'
    
    sio.savemat( os.path.join(out_path, out_fname), {'n_points': observed.shape[0],
                                                     'observed_medaxis_data' :observed_varmean[0], 
                                                     'observed_edge_data'    :observed_varmean[1], 
                                                     'observed_centroid_data':observed_varmean[2],
                                                     'uniform_medaxis_data' :uniform_varmean[0], 
                                                     'uniform_edge_data'    :uniform_varmean[1],
                                                     'uniform_centroid_data':uniform_varmean[2],
                                                     'medaxis_cdf' :[cdf_ma, regions_bd],
                                                     'edge_cdf'    :[cdf_edge, regions_bd],
                                                     'centroid_cdf':[cdf_cent, regions_ct]})



if __name__ == '__main__':

    img_names  = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
                  "blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"]

    conditions = ["bounding_circle","in_shape","touchpoint_hull","patient_fitted"]
    patients   = ["DF","MC"]

    in_path  = "D:/ShapeTapper-Analysis/"
    out_path = "D:/ShapeTapper-Analysis/"

    shapes = ShapeIO(in_path, out_path, img_names)

    shapes.run(spatialAnalysis, patients, conditions)