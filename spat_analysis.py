#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 19 15:42:30 2018

Performs analysis of distance and point clustering (Sadahiro & Takami, 2001) for each shape.

Distance data will be saved in the directory "distance_analysis" with a subdirectory for each shape mask. 
Each file contains data to create a probability distribution defined by the generated uniform sets.

The MATLAB script StatAnalysis.m will use these to calculate the significance of the observed data. 

@author: Jamie Dunkle
"""

import cv2
import os, errno
import numpy as np
import scipy.io as sio
import distances as dist
from timeit import default_timer as timer


################## Globals - CHANGE THESE TO RUN ON SPECIFIC SUBJECTS AND SHAPE SETS
analysis_conds = ["bounding_circle","in_shape","touchpoint_hull"]
img_names = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
             "blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"]
patient = "DF"   
img_path = "./Shapes/"         # path containing shape images
NUM_CDF_STEPS = 100


################# Function Definitions #########################################################################
def load_mat(mat_path, mat_name, img_name) :
    try:
        dat = sio.loadmat(os.path.join(mat_path, mat_name))
    except (TypeError, IOError) :
        print "Error loading: {0} -- skipping {1}...\n".format(mat_name, img_name)
        return None
    return dat


def getVarMean(data) :
    a = data.ndim - 1
    var = np.sum(np.power(data,2),axis=a)/(data.shape[a]-1)
    rmse = np.sqrt(var)
    amd = np.mean(data, axis=a)
    return [var, rmse, amd]


def get_dists(observed, generated, ma_points, edge_points, centroid) :
    #print "Calculating Observed"
    observed_ma_dists   = dist.points2points_c(observed,ma_points)
    observed_edge_dists = dist.points2points_c(observed,edge_points)
    observed_cent_dists = dist.points2points_c(observed,centroid)

    observed_medaxis_data = getVarMean(observed_ma_dists)
    observed_edge_data    = getVarMean(observed_edge_dists)
    observed_cent_data    = getVarMean(observed_cent_dists)
   
    #print "Calculating Generated"
    generated_ma_dists   = np.empty(generated.shape[0:2], dtype=np.float32)
    generated_edge_dists = np.empty(generated.shape[0:2], dtype=np.float32)
    generated_cent_dists = np.empty(generated.shape[0:2], dtype=np.float32)
    
    for i, point_set in enumerate(generated) :
        generated_ma_dists[i]   = dist.points2points_c(point_set,ma_points)
        generated_edge_dists[i] = dist.points2points_c(point_set,edge_points)
        generated_cent_dists[i] = dist.points2points_c(point_set,centroid)

    generated_medaxis_data = getVarMean(generated_ma_dists)
    generated_edge_data    = getVarMean(generated_edge_dists)
    generated_cent_data    = getVarMean(generated_cent_dists)

    observed_dist_data  = [observed_medaxis_data, observed_edge_data, observed_cent_data]
    generated_dist_data = [generated_medaxis_data, generated_edge_data, generated_cent_data]

    return observed_dist_data, generated_dist_data 


def get_cdf(ro_pts, regions, uniform, observed, generated) :
    # calculate d-values for main shape
    ds_uf = dist.cdf_np(uniform, ro_pts, regions)

    # calculate d-values for observed data
    ds_o = dist.cdf_c(observed, ro_pts, regions)

    # calculate d-values for 100k simulated data sets
    ds_g = np.empty((generated.shape[0],regions.shape[0]), dtype=np.float32)
    for i, point_set in enumerate(generated) :
        ds_g[i] = dist.cdf_c(point_set, ro_pts, regions)

    return [ds_uf, ds_o, ds_g]


def matAnalysis(img_name, img_path, img_mat, mat_path, obs_mat, obs_path, gen_mat, gen_path, out_path):    
    ################### Load data ###################

    # Read in the image 
    img_file = img_name + ".png"
    img = cv2.imread(os.path.join(img_path, img_file) ,cv2.IMREAD_UNCHANGED)
    img[(img[:,:,3]==0),0:3] = 0
    
    # Load shape analysis data (medial axis, edge, centroid)
    shape_analysis = load_mat(mat_path, img_mat, img_name)
    if (shape_analysis is None) : return
    ma_points   = np.ascontiguousarray(shape_analysis['ma_points'].astype(np.float32)) # (x,y)
    edge_points = np.ascontiguousarray(shape_analysis['edge_points'].astype(np.float32)) # (x,y)
    centroid    = np.ascontiguousarray(shape_analysis['centroid'].astype(np.float32)) # (x,y)
    
    # Load observed touchpoint data
    observed_mat = load_mat(obs_path, obs_mat, img_name)
    if (observed_mat is None) : return
    observed = np.ascontiguousarray(observed_mat['img_dataset'].astype(np.float32))
    observed[:,1] = img.shape[0]-observed[:,1] # opencv coordinates use inverted y-axis (top-left origin)
    
    # Load generated uniform data
    generated_mat = load_mat(gen_path, gen_mat, img_name)
    if (generated_mat is None) : return
    generated = np.ascontiguousarray(generated_mat['gen_datasets'].astype(np.float32))

    # Trim observed touchpoints outside of shape for the "in shape" condition
    if cond == "in_shape" :
        observed_inshape = [] 
        for i in range(observed.shape[0]) :
            if cv2.pointPolygonTest(edge_points, tuple(observed[i]), measureDist=False) == 1 :
                observed_inshape.append(observed[i])
        observed = np.array(observed_inshape)

    ################### Spatial Analysis ################### 

    # Calculate variance & mean info for observed and generated sets
    start = timer()
    observed_varmean, generated_varmean = get_dists(observed, generated, ma_points, edge_points, centroid)
    print "var-mean data :", timer() - start, "s"
    
    # Find uniformly distributed points in shape (pixel locations)
    img_bin = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
    img_bin[(img_bin!=0)] = 1 # convert to logical array (white = 1)
    white_idx = np.where(img_bin==1)    
    shape_points = np.flip(np.transpose(np.array(white_idx)), axis=1)
    shape_points = np.ascontiguousarray(shape_points, dtype=np.float32)

    # Calculate vector of regions for spatial analysis (from 0 to bounds of enclosing circle)
    # Number of steps is specified in NUM_RS variable at top of script
    img_bounds = np.array( [ [0,0], [img.shape[0], 0] , [img.shape[0],img.shape[1]], [0,img.shape[1]] ])
    _, bd_circ_radius = cv2.minEnclosingCircle(img_bounds)
    max_r_bd = int(np.ceil(bd_circ_radius))
    regions_bd = np.linspace(0.0, max_r_bd, num=NUM_CDF_STEPS, endpoint=True)
    regions_bd = np.ascontiguousarray(regions_bd, dtype=np.float32)
    
    max_r_ct = int(np.ceil(bd_circ_radius*np.sqrt(2)))
    regions_ct = np.linspace(0.0, max_r_ct, num=NUM_CDF_STEPS, endpoint=True)
    regions_ct = np.ascontiguousarray(regions_ct, dtype=np.float32)

    # Calculate spatial CDF data for each Reference Object (Medax, Edge, Centroid)
    start = timer()
    cdf_ma   = get_cdf(ma_points, regions_ct, shape_points, observed, generated)
    cdf_edge = get_cdf(edge_points, regions_bd, shape_points, observed, generated)
    cdf_cent = get_cdf(centroid, regions_ct, shape_points, observed, generated)
    print "spatial data  :", timer()-start, "s\n"
    
    # Output data for statistical analysis
    sio.savemat(os.path.join(out_path,img_name+'_analysis.mat'),
                {'n_points': observed.shape[0],
                 'observed_medaxis_data' :observed_varmean[0], 
                 'observed_edge_data'    :observed_varmean[1], 
                 'observed_centroid_data':observed_varmean[2],
                 'generated_medaxis_data' :generated_varmean[0], 
                 'generated_edge_data'    :generated_varmean[1],
                 'generated_centroid_data':generated_varmean[2],
                 'medaxis_cdf' :[cdf_ma, regions_ct],
                 'edge_cdf'    :[cdf_edge, regions_bd],
                 'centroid_cdf':[cdf_cent, regions_ct]})
    

if __name__ == '__main__':
    
    mat_path = img_path+"shape_analysis/"          # path containing medial axis mat files
    obs_path = patient+"/aggregated_observations/" # path containing observed data
    
    for cond in analysis_conds :
        print "Condition:", cond, '\n'
        
        gen_path = patient+"/generated_uniform_data/"+cond+"/"  # path containing generated uniform data
        out_path = patient+"/distance_analysis/"+cond+"/"       # output path

        try:
            os.makedirs(out_path)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

        for img_name in img_names :
            print 'Starting', img_name
            img_mat = img_name + "_shape_analysis.mat"
            obs_mat = img_name + "_Patient_"+patient+"_aggregated_observations.mat"
            gen_mat = img_name + "_Patient_"+patient+"_generated_uniform_sets.mat"
            matAnalysis(img_name, img_path, img_mat, mat_path, obs_mat, obs_path, gen_mat, gen_path, out_path)