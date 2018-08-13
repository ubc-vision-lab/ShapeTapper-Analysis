#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Performs analysis of distance and point clustering (Sadahiro & Takami, 2001) for each shape.

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


# Gets Average Mean Distance and StdDev for points to each reference object
def get_dists(points, medial_axis, edge_points, centroid) :
    if points is None or len(points.shape) <= 1 : return [[],[],[]]

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


# Cumulative Density Function for "Reference Object" defined by ro_pts
def cdf(points_in, ro_pts, regions):
    ft_dists = np.sort( dist.points2points(points_in,ro_pts) ) # sorted distances from RO
    cdf = np.empty(regions.shape[0], dtype=np.float32) # output vector
    gu_cdf(regions, ft_dists, cdf) # count number of points in each region
    return np.true_divide(cdf, points_in.shape[0]) # normalize [0, n_pts] to [0,1]


# Gets CDFs for reference object (interior (uniform pixel) points = expected)
def get_cdf(ro_pts, regions, interior_points, observed, uniform=None) :
    # calculate CDF for main shape
    ds_exp = cdf(interior_points, ro_pts, regions)

    # calculate CDF for observed data, get difference from expected
    if len(observed.shape) > 1 :
        ds_obs = cdf(observed, ro_pts, regions)
        obs_diff = ds_obs - ds_exp
    else :
        obs_diff = [0]
    
    obs_dplus = np.max(obs_diff)
    obs_dminus = np.min(obs_diff)
    obs_dplus_r = np.argmax(obs_diff)
    obs_dminus_r = np.argmin(obs_diff)
    
    # calculate CDF for 100k simulated data sets, get difference from expected
    if uniform is not None :
        ds_unif = np.empty((uniform.shape[0],regions.shape[0]), dtype=np.float32)
        for i, point_set in enumerate(uniform) :
            ds_unif[i] = cdf(point_set, ro_pts, regions)
        gen_diff = ds_unif - ds_exp
        gen_dplus = np.max(gen_diff, axis=1)
        gen_dminus = np.min(gen_diff, axis=1)
    else :
        gen_dplus = [0]
        gen_dminus = [0]
    
    return [obs_dplus, obs_dminus, obs_dplus_r, obs_dminus_r, gen_dplus, gen_dminus]



def spatialAnalysis(shape, out_path, patient, cond, task=None, n_cdf=1000):    

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

    if task is None :
        out_path = os.path.join(out_path, patient, "spatial_analysis", cond)
    else :
        out_path = os.path.join(out_path, patient, "spatial_analysis", cond, task)

    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    # Fetch observed and generated uniform touch points
    observed = shape.observed
    uniform  = shape.uniform

    # Fetch Reference Object points
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
    # Number of steps defaults to 1000
    img_bounds = np.array( [ [0,0], [shape.dims[0], 0] , [shape.dims[0],shape.dims[1]], [0,shape.dims[1]] ])
    _, bd_circ_radius = cv2.minEnclosingCircle(img_bounds)
    max_r_bd = int(np.ceil(bd_circ_radius))
    regions_bd = np.linspace(0.0, max_r_bd, num=n_cdf, endpoint=True)
    regions_bd = np.ascontiguousarray(regions_bd, dtype=np.float32)
    
    # Use a slightly enlarged regions vector for centroid to catch all points
    max_r_ct = int(np.ceil(bd_circ_radius*np.sqrt(2))) 
    regions_ct = np.linspace(0.0, max_r_ct, num=n_cdf, endpoint=True)
    regions_ct = np.ascontiguousarray(regions_ct, dtype=np.float32)

    # Calculate spatial CDF data for each Reference Object (Medaxis, Edge, Centroid)
    start = timer()
    cdf_ma   = get_cdf(medial_axis, regions_bd, shape_points, observed, uniform)
    cdf_edge = get_cdf(edge_points, regions_bd, shape_points, observed, uniform)
    cdf_cent = get_cdf(centroid, regions_ct, shape_points, observed, uniform)
    print "{0} analysis : cdf data generated in {1}s".format(shape.name, timer() - start)

    # Generate MAT file output names
    if shape.pair_mapping is not None :
        if task is None :
            out_fname = "_".join((shape.pair_mapping, "to", shape.name, "Patient", patient, "spatial_analysis", cond)) + '.mat'
        else :
            out_fname = "_".join((shape.pair_mapping, "to", shape.name, "Patient", patient, "spatial_analysis", cond, task)) + '.mat'
    else :
        if task is None :
            out_fname = "_".join((shape.name, "Patient", patient, "spatial_analysis", cond)) + '.mat'
        else :
            out_fname = "_".join((shape.name, "Patient", patient, "spatial_analysis", cond, task)) + '.mat'

    # Output data for statistical analysis in the StatAnalysis.m scripts
    sio.savemat( os.path.join(out_path, out_fname), {'n_points'              : observed.shape[0],
                                                     'observed_medaxis_data' : observed_varmean[0], 
                                                     'observed_edge_data'    : observed_varmean[1], 
                                                     'observed_centroid_data': observed_varmean[2],
                                                     'uniform_medaxis_data'  : uniform_varmean[0], 
                                                     'uniform_edge_data'     : uniform_varmean[1],
                                                     'uniform_centroid_data' : uniform_varmean[2],
                                                     'medaxis_cdf_obs' : [cdf_ma[0], cdf_ma[1]],
                                                     'medaxis_cdf_obs_r' : [cdf_ma[2], cdf_ma[3]],
                                                     'medaxis_cdf_gen' : [cdf_ma[4], cdf_ma[5]],
                                                     'edge_cdf_obs' : [cdf_edge[0], cdf_edge[1]],
                                                     'edge_cdf_obs_r' : [cdf_edge[2], cdf_edge[3]],
                                                     'edge_cdf_gen' : [cdf_edge[4], cdf_edge[5]],
                                                     'centroid_cdf_obs' : [cdf_cent[0], cdf_cent[1]],
                                                     'centroid_cdf_obs_r' : [cdf_cent[2], cdf_cent[3]],
                                                     'centroid_cdf_gen' : [cdf_cent[4], cdf_cent[5]],
                                                     'regions_bd' : regions_bd,
                                                     'regions_ct' : regions_ct })



if __name__ == '__main__':

    img_names  = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
                  "blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"]

    conditions = ["bounding_circle","in_shape","touchpoint_hull","patient_fitted"]
    patients   = ["DF","MC"]

    in_path  = "D:/ShapeTapper-Analysis/"
    out_path = "D:/ShapeTapper-Analysis/"

    shapes = ShapeIO(in_path, out_path, img_names)

    shapes.run(spatialAnalysis, patients, conditions)