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
             "blake_01","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"]
patient = "MC"   
img_path = "./Shapes/"         # path containing shape images



################# Function Definitions #########################################################################
def getVarMean(data) :
    var = np.sum(np.power(data, 2))/(np.size(data)-1)
    rmse = np.sqrt(var)
    amd = np.mean(data)
    return [var, rmse, amd]


def getVarMean2D(data) :
    var = np.sum(np.power(data,2),axis=1) / (data.shape[1] - 1)
    rmse = np.sqrt(var)
    amd = np.mean(data, axis=1)
    return [var, rmse, amd]


# # cumulative distribution function for "reference object" defined by ro_pts
# # returns ratio of points located within a given distance of the reference object
# def cdf(points_in, ro_pts, max_r):
#     # points2points_cmin crashes for points > 250, so fall back to numpy function
#     if points_in.shape[0] > 250 : 
#         ft_dists = dist.points2points_np(points_in,ro_pts)
#     else :
#         ft_dists = dist.points2points_c(points_in,ro_pts)
#     points_in_region = np.empty((max_r+1,1), dtype=int)
#     for i in range(max_r+1) : 
#         points_in_region[i] = (ft_dists <= i).sum()
#     return np.true_divide(points_in_region, points_in.shape[0])


def get_dvals(ro_pts, max_r, generated_unif, observed, generated) :
    # calculate d-values for main shape
    ds_gu = dist.cdf(generated_unif, ro_pts, max_r)
    # calculate d-values for observed data
    ds_o = dist.cdf(observed, ro_pts, max_r)
    diff_o = ds_o - ds_gu
    d_o_plus = np.max(diff_o)
    d_o_plus_r = np.argmax(diff_o)
    d_o_minus = np.min(diff_o)
    d_o_minus_r = np.argmin(diff_o)
    # calculate d-values for 100k simulated data sets --> CDF
    ds_g = np.zeros((generated.shape[0],max_r), dtype=np.float32)
    diff_g = np.zeros((generated.shape[0],max_r), dtype=np.float32)
    for i in range(generated.shape[0]) :
        ds_g[i,:] = dist.cdf(generated[i], ro_pts, max_r)
        diff_g[i,:] = ds_g[i,:] - ds_gu
    d_g_plus = np.max(diff_g, axis=1)
    d_g_minus = np.min(diff_g, axis=1)
    return d_o_plus, d_o_plus_r, d_o_minus, d_o_minus_r, d_g_plus, d_g_minus


def matAnalysis(img_file, img_path, img_mat, mat_path, obs_mat, obs_path, gen_mat, gen_path, out_path):    

    # Read in the image 
    img = cv2.imread(os.path.join(img_path, img_file) ,cv2.IMREAD_UNCHANGED)
    img[(img[:,:,3]==0),0:3] = 0
    
    try:
        shape_analysis = sio.loadmat(os.path.join(mat_path, img_mat))
    except (TypeError, IOError) :
        print "Error loading: {0} -- skipping {1}...".format(img_mat, img_name)
        return
 
    ma_points = np.ascontiguousarray(shape_analysis['ma_points'].astype(np.float32)) # (x,y)
    edge_points = np.ascontiguousarray(shape_analysis['edge_points'].astype(np.float32)) # (x,y)
    centroid = np.ascontiguousarray(shape_analysis['centroid'].astype(np.float32)) # (x,y)
    
    #print "Calculating Observed"
    try:
        observed_mat = sio.loadmat(os.path.join(obs_path, obs_mat))
    except (TypeError, IOError) :
        print "Error loading: {0} -- skipping {1}...".format(obs_mat, img_name)
        return
    observed = np.ascontiguousarray(observed_mat['img_dataset'].astype(np.float32))
    observed[:,1] = img.shape[0]-observed[:,1] # opencv coordinates use inverted y-axis
    
    if cond == "in_shape" :
        observed_inshape = [] 
        for i in range(observed.shape[0]) :
            if cv2.pointPolygonTest(edge_points, tuple(observed[i]), measureDist=False) == 1 :
                observed_inshape.append(observed[i])
        observed = np.array(observed_inshape)
    
    start = timer()
    observed_ma_dists = dist.points2points_c(observed,ma_points)
    observed_edge_dists = dist.points2points_c(observed,edge_points)
    observed_centroid_dists = dist.points2points_c(observed,centroid)

     # Save data
    observed_medaxis_data = getVarMean(observed_ma_dists)
    observed_edge_data = getVarMean(observed_edge_dists)
    observed_centroid_data = getVarMean(observed_centroid_dists)
   
#    #print "Calculating Generated"
    generated_mat = sio.loadmat(os.path.join(gen_path, gen_mat))
    generated = np.ascontiguousarray(generated_mat['gen_datasets'].astype(np.float32))
    
    generated_ma_dists = np.empty((generated.shape[0],generated.shape[1]))
    generated_edge_dists = np.empty((generated.shape[0],generated.shape[1]))
    generated_centroid_dists = np.empty((generated.shape[0],generated.shape[1]))
    
    for i, point_set in enumerate(generated) :
        generated_ma_dists[i] = dist.points2points_c(point_set,ma_points)
        generated_edge_dists[i] = dist.points2points_c(point_set,edge_points)
        generated_centroid_dists[i] = dist.points2points_c(point_set,centroid)

    #Save data
    generated_medaxis_data = getVarMean2D(generated_ma_dists)
    generated_edge_data = getVarMean2D(generated_edge_dists)
    generated_centroid_data = getVarMean2D(generated_centroid_dists)
    generated_medaxis_data = [list(datum) for datum in generated_medaxis_data]
    generated_edge_data = [list(datum) for datum in generated_edge_data]
    generated_centroid_data = [list(datum) for datum in generated_centroid_data]
    print "var-mean data:", timer() - start, "s"
    
    sio.savemat(os.path.join(out_path,img_name+'_analysis_var_mean.mat'),
                {'observed_medaxis_data':observed_medaxis_data, 
                 'observed_edge_data':observed_edge_data, 
                 'observed_centroid_data':observed_centroid_data,
                 'generated_medaxis_data':generated_medaxis_data, 
                 'generated_edge_data':generated_edge_data,
                 'generated_centroid_data':generated_centroid_data})
    
    img_bin = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
    img_bin[(img_bin!=0)] = 1 # convert to logical array (white = 1)
    white_idx = np.where(img_bin==1)    
    shape_points = np.transpose(np.array(white_idx))
    shape_points = np.ascontiguousarray(np.flip(shape_points, axis=1).astype(np.float32))
    
    max_r = np.max(img.shape[:2])/2
    
    start = timer()
    dm_o_plus, dm_o_plus_r, dm_o_minus, dm_o_minus_r, dm_g_plus, dm_g_minus = get_dvals(ma_points, max_r, shape_points, observed, generated)
    print "spat medial axis data:", timer() - start, "s"
    
    start = timer()
    de_o_plus, de_o_plus_r, de_o_minus, de_o_minus_r, de_g_plus, de_g_minus = get_dvals(edge_points, max_r, shape_points, observed, generated)
    print "spat edge data:", timer() - start, "s"
    
    start = timer()
    dc_o_plus, dc_o_plus_r, dc_o_minus, dc_o_minus_r, dc_g_plus, dc_g_minus = get_dvals(centroid, max_r, shape_points, observed, generated)
    print "spat centroid data:", timer()-start, "s\n"
    
    sio.savemat(os.path.join(out_path,img_name+'_analysis_spat.mat'),
                {'observed_medaxis_dplus':[dm_o_plus,dm_o_plus_r],
                 'observed_medaxis_dminus':[dm_o_minus,dm_o_minus_r],
                 'observed_edge_dplus':[de_o_plus,de_o_plus_r],
                 'observed_edge_dminus':[de_o_minus,de_o_minus_r],
                 'observed_centroid_dplus':[dc_o_plus,dc_o_plus_r],
                 'observed_centroid_dminus':[dc_o_minus,dc_o_minus_r],
                 'generated_medaxis_dplus':dm_g_plus,
                 'generated_medaxis_dminus':dm_g_minus,
                 'generated_edge_dplus':de_g_plus,
                 'generated_edge_dminus':de_g_minus,
                 'generated_centroid_dplus':dc_g_plus,
                 'generated_centroid_dminus':dc_g_minus})
    
    

if __name__ == '__main__':
    
    mat_path = patient+"/shape_analysis/"          # path containing medial axis mat files
    obs_path = patient+"/aggregated_observations/" # path containing observed data
    
    for cond in analysis_conds :
        print "Condition:", cond
        
        gen_path = patient+"/generated_uniform_data/"+cond+"/"  # path containing generated uniform data
        
        out_path = patient+"/distance_analysis/"+cond+"/"       # output path
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
            gen_mat = img_name + "_Patient_"+patient+"_generated_uniform_sets.mat"
            matAnalysis(img_file, img_path, img_mat, mat_path, obs_mat, obs_path, gen_mat, gen_path, out_path)