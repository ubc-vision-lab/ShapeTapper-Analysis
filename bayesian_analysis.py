#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Prototype Bayesian analysis comparing the likelihood of the Reference Objects to explain observed touch points

@author: Jamie Dunkle
"""

import cv2
import sys
import os, errno
import numpy as np
import scipy.io as sio
from scipy.ndimage.filters import gaussian_filter

################# Function Definitions #########################################################################

def make_gauss_prob(points, dims, sigma) :
    pt_array = np.zeros(dims, dtype=np.float32)
    pt_array[zip(*points)] = 1
    pt_filt = gaussian_filter(pt_array, sigma, mode='constant', cval=0.0)
    return pt_filt / np.sum(pt_filt)


def bayesAnalysis(shape, out_path, patient, cond=None, max_sigma = 100):    

    img = shape.img
    observed = shape.observed

    medial_axis = shape.medial_axis
    edge_points = shape.edge_points
    centroid    = shape.centroid

    bounds = np.array([ [0,0], [shape.dims[1],0], [shape.dims[1],shape.dims[0]], [0,shape.dims[0]] ])
    
    # Find all touchpoints in image canvas
    observed_incanvas = []
    oob = 0
    for i in range(observed.shape[0]) :
        if cv2.pointPolygonTest(bounds, tuple(observed[i]), measureDist=False) == 1 :
            observed_incanvas.append(observed[i])
        else :
            oob += 1
    observed = np.array(observed_incanvas)
    
    # Create null hypothesis: a normalized canvas with identical uniform values
    null_prob  = np.ones(img.shape[1::-1], dtype=np.float32)
    null_prob /= np.sum(null_prob)

    # Find prob of null, sum of null hypothesis at all observed points
    p_null = np.sum(null_prob[zip(*observed)])
    
    # Calculate hypotheses for each Reference Object (experimental sigma values)
    p_mapt = []
    p_edge = []
    p_cent = []
    for i in range(1,max_sigma+1):
        if i%10 == 0 : sys.stdout.write(".")
        mapt_filt = make_gauss_prob(medial_axis, img.shape[1::-1], sigma=i)
        edge_filt = make_gauss_prob(edge_points, img.shape[1::-1], sigma=i)
        cent_filt = make_gauss_prob(centroid, img.shape[1::-1], sigma=i)        
        p_mapt.append(np.sum(mapt_filt[zip(*observed)]))
        p_edge.append(np.sum(edge_filt[zip(*observed)]))
        p_cent.append(np.sum(cent_filt[zip(*observed)]))
    
    # Find sigma value for which prob of all hypothesis is greatest
    max_idx = np.argmax(np.sum(np.array([p_cent, p_edge, p_mapt]),0))
    
    # Print calculated probs at max sigma
    print "\nNull hypothesis :", p_null
    print "Centroid :", p_cent[max_idx]
    print "Edge :", p_edge[max_idx]
    print "MedAx :", p_mapt[max_idx], '\n'
    
    # Output data for all sigmas to MAT file
    out_fname = os.path.join(out_path, shape.name+'_bayes.mat')
    sio.savemat( out_fname ,
                {'null_p': p_null,
                 'centroid_p':[p_cent, range(1,max_sigma+1)], 
                 'edge_p':[p_edge, range(1,max_sigma+1)], 
                 'medaxis_p':[p_mapt, range(1,max_sigma+1)]})
    

if __name__ == '__main__':
    pass