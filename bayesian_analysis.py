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
import sys
import os, errno
import numpy as np
import scipy.io as sio
from scipy.ndimage.filters import gaussian_filter


################## Globals - CHANGE THESE TO RUN ON SPECIFIC SUBJECTS AND SHAPE SETS
img_names = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
             "blake_01","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"]
patient = "MC"   
img_path = "./Shapes/"         # path containing shape images



################# Function Definitions #########################################################################

def make_gauss_prob(points, dims, sigma) :
    pt_array = np.zeros(dims, dtype=np.float32)
    pt_array[zip(*points)] = 1
    pt_filt = gaussian_filter(pt_array, sigma, mode='constant', cval=0.0)
    return pt_filt / np.sum(pt_filt)


def matAnalysis(img_file, img_path, img_mat, mat_path, obs_mat, obs_path, out_path):    

    # Read in the image 
    img = cv2.imread(os.path.join(img_path, img_file) ,cv2.IMREAD_UNCHANGED)
    img[(img[:,:,3]==0),0:3] = 0
    
    try:
        shape_analysis = sio.loadmat(os.path.join(mat_path, img_mat))
    except (TypeError, IOError) :
        print "Error loading: {0} -- skipping {1}...".format(img_mat, img_name)
        return
 
    ma_points = shape_analysis['ma_points'].astype(int) # (x,y)
    edge_points = shape_analysis['edge_points'].astype(int) # (x,y)
    centroid = shape_analysis['centroid'].astype(int) # (x,y)
    
    #print "Calculating Observed"
    try:
        observed_mat = sio.loadmat(os.path.join(obs_path, obs_mat))
    except (TypeError, IOError) :
        print "Error loading: {0} -- skipping {1}...".format(obs_mat, img_name)
        return
    observed = np.ascontiguousarray(observed_mat['img_dataset'].astype(int))
    observed[:,1] = img.shape[0]-observed[:,1] # opencv coordinates use inverted y-axis

    bounds = np.array([[0,0],[img.shape[1],0],[img.shape[1],img.shape[0]],[0,img.shape[0]]])
    observed_incanvas = []
    oob = 0
    for i in range(observed.shape[0]) :
        if cv2.pointPolygonTest(bounds, tuple(observed[i]), measureDist=False) == 1 :
            observed_incanvas.append(observed[i])
        else :
            oob += 1
    observed = np.array(observed_incanvas)
    
    null_prob  = np.ones(img.shape[1::-1], dtype=np.float32)
    null_prob /= np.sum(null_prob)
    p_null = np.sum(null_prob[zip(*observed)])
    
    p_mapt = []
    p_edge = []
    p_cent = []
    for i in range(1,101):
        if i%10 == 0 : sys.stdout.write(".")
        mapt_filt = make_gauss_prob(ma_points, img.shape[1::-1], sigma=i)
        edge_filt = make_gauss_prob(edge_points, img.shape[1::-1], sigma=i)
        cent_filt = make_gauss_prob(centroid, img.shape[1::-1], sigma=i)        
        p_mapt.append(np.sum(mapt_filt[zip(*observed)]))
        p_edge.append(np.sum(edge_filt[zip(*observed)]))
        p_cent.append(np.sum(cent_filt[zip(*observed)]))
        
    max_idx = np.argmax(np.sum(np.array([p_cent, p_edge, p_mapt]),0))
    
    print "\nNull hypothesis :", p_null
    print "Centroid :", p_cent[max_idx]
    print "Edge :", p_edge[max_idx]
    print "MedAx :", p_mapt[max_idx], '\n'
    
    sio.savemat(os.path.join(out_path,img_name+'_bayes.mat'),
                {'null_p': p_null,
                 'centroid_p':[p_cent, range(1,101)], 
                 'edge_p':[p_edge, range(1,101)], 
                 'medaxis_p':[p_mapt, range(1,101)]})
    

    
    

if __name__ == '__main__':
    
    mat_path = img_path+"shape_analysis/"          # path containing medial axis mat files
    obs_path = patient+"/aggregated_observations/" # path containing observed data
        
    out_path = patient+"/distance_analysis/"       # output path
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
        matAnalysis(img_file, img_path, img_mat, mat_path, obs_mat, obs_path, out_path)