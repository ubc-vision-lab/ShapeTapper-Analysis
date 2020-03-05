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


def bayesianAnalysis(shape, out_path, patient, cond=None, task=None, max_sigma=80):    
    
    if patient is None :
        print("Error in bayesianAnalysis() : no patient name specified.")
        return

    if shape.observed is None :
        print("Error in bayesianAnalysis() : no observed touch points found for {0}. Please check data files.".format(shape.name))
        return

    out_path = os.path.join(out_path, patient, "bayesian_analysis")
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    print("Starting Bayes analysis on {0}".format(shape.name))

    img = np.copy(shape.img)
    observed = shape.observed.astype(np.int32)

    medial_axis = shape.medial_axis.astype(np.int32)
    edge_points = shape.edge_points.astype(np.int32)
    centroid    = shape.centroid.astype(np.int32)

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
    # null_prob  = np.ones(img.shape[1::-1], dtype=np.float32)
    # null_prob /= np.sum(null_prob)
    null_prob = np.copy(img[:,:,3])
    null_prob = np.transpose(null_prob)
    null_prob[null_prob==255] = 1
    null_prob = null_prob.astype(np.float32)
    null_prob /= np.sum(null_prob)

    # Find prob of null, sum of null hypothesis at all observed points
    p_null = np.sum(null_prob[zip(*observed)])
    
    # Calculate hypotheses for each Reference Object (experimental sigma values)
    p_mapt = []
    p_mapt_zval = []
    # p_mapt_sigma = []
    # p_edge = []
    p_cent = []
    p_cent_zval = []
    p_cent_sigma = []
    p_mapt_plus_cent = []
    #p_mapt_solo = []
    #p_cent_solo = []
    for i in range(0,max_sigma):
# =============================================================================
#         if i%10 == 0 : sys.stdout.write(".")
#         # mapt_filt = make_gauss_prob(medial_axis, img.shape[1::-1], sigma=i)
#         mapt_filt = make_gauss_prob(medial_axis, img.shape[1::-1], sigma=i)
#         mapt_filt_zval = np.max(mapt_filt[zip(*medial_axis)])
#         found_ct_sigma = False
#         ct_sigma = 0.0
#         while not found_ct_sigma :
#             cent_filt = make_gauss_prob(centroid, img.shape[1::-1], sigma=ct_sigma)
#             # sys.stdout.write(".")
#             cent_filt_zval = np.max(cent_filt[zip(*centroid)])
#             if mapt_filt_zval > cent_filt_zval :
#                 found_ct_sigma = True
#             else :
#                 ct_sigma += 0.05
# =============================================================================
        
        in_path = 'C:\\ShapeTapper-Analysis\\Shapes\\bayesian_shape_analysis\\'
        in_file = "_".join((shape.name, "gauss_ma_sigma", str(i))) + '.mat'
        full_in_path = os.path.join(in_path, in_file)
        try:
            shape_gauss_mat = sio.loadmat(full_in_path)
        except (TypeError, IOError) :
            print "Error loading bayes shape data: {0} -- skipping...".format(full_in_path)
            return
        
        mapt_filt = np.ascontiguousarray(shape_gauss_mat['mapt_filt'], dtype=np.float32)
        mapt_filt_zval = np.ascontiguousarray(shape_gauss_mat['mapt_filt_zval'], dtype=np.float32)
        cent_filt = np.ascontiguousarray(shape_gauss_mat['cent_filt'], dtype=np.float32)
        cent_filt_zval = np.ascontiguousarray(shape_gauss_mat['cent_filt_zval'], dtype=np.float32)
        ct_sigma = np.ascontiguousarray(shape_gauss_mat['sigma_ct'], dtype=np.float32)
        
# =============================================================================
#         out_path = os.path.join("C:\\ShapeTapper-Analysis\\", "Shapes", "bayesian_shape_analysis")
#         try:
#             os.makedirs(out_path)
#         except OSError as e:
#             if e.errno != errno.EEXIST:
#                 raise
#         out_fname = "_".join((shape.name, "gauss_ma_sigma", str(i))) + '.mat'
#         
#         sio.savemat( os.path.join(out_path, out_fname) ,
#             {'sigma_ma'    : i,
#              'sigma_ct': ct_sigma, 
#              'mapt_filt': mapt_filt,
#              'mapt_filt_zval' : mapt_filt_zval,
#              'cent_filt' : cent_filt,
#              'cent_filt_zval' : cent_filt_zval} )
# =============================================================================
    
        # edge_filt = make_gauss_prob(edge_points, img.shape[1::-1], sigma=i)
        # cent_filt = make_gauss_prob(centroid, img.shape[1::-1], sigma=i)
        p_mapt.append(np.sum(mapt_filt[zip(*observed)]))
        p_mapt_zval.append(mapt_filt_zval)
        #p_mapt_sigma.append(ct_sigma)
        # p_edge.append(np.sum(edge_filt[zip(*observed)]))
        p_cent.append(np.sum(cent_filt[zip(*observed)]))
        p_cent_zval.append(cent_filt_zval)
        p_cent_sigma.append(ct_sigma)
        
        mapt_plus_cent_filt = mapt_filt + cent_filt
        mapt_plus_cent_filt /= np.sum(mapt_plus_cent_filt)
        p_mapt_plus_cent.append(np.sum(mapt_plus_cent_filt[zip(*observed)]))
        #mapt_solo_filt = mapt_filt - cent_filt
        #mapt_solo_filt[mapt_solo_filt < 0] = 0
        #cent_solo_filt = cent_filt - mapt_solo_filt
        #p_cent_solo.append(np.sum(cent_solo_filt[zip(*observed)]))
        #mapt_sum = np.sum(mapt_solo_filt) 
        #if mapt_sum != 0 :
        #    mapt_solo_filt /= mapt_sum      
        #else :
        #    mapt_solo_filt  = np.zeros_like(mapt_solo_filt)
        #p_mapt_solo.append(np.sum(mapt_solo_filt[zip(*observed)]))
        
    print('\n')

    # Find sigma value for which prob of all hypothesis is greatest
    # max_idx = np.argmax(np.sum(np.array([p_cent, p_edge, p_mapt]),0))
    max_idx = np.argmax(np.sum(np.array([p_cent, p_mapt]),0))
    
    # Print calculated probs at max sigma
    print("Null hypothesis : {0}".format(p_null))
    print("Centroid : {0}".format(p_cent[max_idx]))
    #print("Edge : {0}".format(p_edge[max_idx]))
    print("MedAx : {0}\n".format(p_mapt[max_idx]))
    
    # Output data for all sigmas to MAT file
    if task is None:
        out_fname = "_".join((shape.name, "Patient", patient, "bayes_rev")) + '.mat'
    else:
        out_fname = "_".join((shape.name, "Patient", patient, task, "bayes_rev")) + '.mat'
    
    #sio.savemat( os.path.join(out_path, out_fname) ,
    sio.savemat( os.path.join('C:\\ShapeTapper-Analysis\\bayes_all_parts\\', out_fname) ,
                {'null_p'    : p_null,
                 'centroid_p': p_cent, 
                 'centroid_zmax': p_cent_zval,
                 'centroid_sigma' : p_cent_sigma,
                 'medaxis_p' : p_mapt,
                 'medaxis_zmax' : p_mapt_zval,
                 'ma_plus_cent_p' : p_mapt_plus_cent} )
                 #'medaxis_sigma' : p_mapt_sigma,
                 #'medaxis_solo_p' : p_mapt_solo 
                 #'centroid_solo_p': p_cent_solo,
                 #'edge_p'    : p_edge 


def bayesianAnalysisUnifDiffScores(shape, out_path, patient, cond=None, task=None, sigma=8):    
    
    if patient is None :
        print("Error in bayesianAnalysis() : no patient name specified.")
        return

    if shape.observed is None :
        print("Error in bayesianAnalysis() : no observed touch points found for {0}. Please check data files.".format(shape.name))
        return

    out_path = os.path.join(out_path, patient, "bayesian_analysis")
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    print("Starting Bayes analysis on {0}".format(shape.name))

    img = np.copy(shape.img)
    # Fetch observed and generated uniform touch points
    observed = shape.observed.astype(np.int32)
    uniform  = shape.uniform.astype(np.int32)

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
    # null_prob  = np.ones(img.shape[1::-1], dtype=np.float32)
    # null_prob /= np.sum(null_prob)
    null_prob = np.copy(img[:,:,3])
    null_prob = np.transpose(null_prob)
    null_prob[null_prob==255] = 1
    null_prob = null_prob.astype(np.float32)
    null_prob /= np.sum(null_prob)

    # Find prob of null, sum of null hypothesis at all observed points
    p_null = np.sum(null_prob[zip(*observed)])
    
    # Load hypotheses for each Reference Object (defined sigma value)
    in_path = 'C:\\ShapeTapper-Analysis\\Shapes\\bayesian_shape_analysis\\'
    in_file = "_".join((shape.name, "gauss_ma_sigma", str(sigma))) + '.mat'
    full_in_path = os.path.join(in_path, in_file)
    try:
        shape_gauss_mat = sio.loadmat(full_in_path)
    except (TypeError, IOError) :
        print "Error loading bayes shape data: {0} -- skipping...".format(full_in_path)
        return
    
    mapt_filt = np.ascontiguousarray(shape_gauss_mat['mapt_filt'], dtype=np.float32)
    mapt_filt_zval = np.ascontiguousarray(shape_gauss_mat['mapt_filt_zval'], dtype=np.float32)
    cent_filt = np.ascontiguousarray(shape_gauss_mat['cent_filt'], dtype=np.float32)
    cent_filt_zval = np.ascontiguousarray(shape_gauss_mat['cent_filt_zval'], dtype=np.float32)
    ct_sigma = np.ascontiguousarray(shape_gauss_mat['sigma_ct'], dtype=np.float32)
    
    # edge_filt = make_gauss_prob(edge_points, img.shape[1::-1], sigma=i)
    # cent_filt = make_gauss_prob(centroid, img.shape[1::-1], sigma=i)
    p_mapt_obs = np.sum(mapt_filt[zip(*observed)])
    p_mapt_zval = mapt_filt_zval
    #p_mapt_sigma.append(ct_sigma)
    # p_edge.append(np.sum(edge_filt[zip(*observed)]))
    p_cent_obs = np.sum(cent_filt[zip(*observed)])
    p_cent_zval = cent_filt_zval
    
    cent_sigma = ct_sigma
    
    p_mapt_unif = np.empty(uniform.shape[0], dtype=np.float32)
    p_cent_unif = np.empty_like(p_mapt_unif)
    for i, point_set in enumerate(uniform) :
        p_mapt_unif[i] = np.sum(mapt_filt[zip(*point_set)])
        p_cent_unif[i] = np.sum(cent_filt[zip(*point_set)])
        
    
    # Output data for all sigmas to MAT file
    if task is None:
        out_fname = "_".join((shape.name, "Patient", patient, "bayes_unif_diff")) + '.mat'
    else:
        out_fname = "_".join((shape.name, "Patient", patient, task, "bayes_unif_diff")) + '.mat'
    
    #sio.savemat( os.path.join(out_path, out_fname) ,
    #sio.savemat( os.path.join('C:\\ShapeTapper-Analysis\\bayes_all_parts\\', out_fname) ,
    sio.savemat( os.path.join(out_path, out_fname) ,
                {'null_p'    : p_null,
                 'centroid_p_obs': p_cent_obs, 
                 'centroid_zmax': p_cent_zval,
                 'centroid_sigma' : cent_sigma,
                 'medaxis_p_obs' : p_mapt_obs,
                 'medaxis_zmax' : p_mapt_zval,
                 'centroid_p_unif': p_cent_unif,
                 'medaxis_p_unif': p_mapt_unif } )
                 #'medaxis_sigma' : p_mapt_sigma,
                 #'medaxis_solo_p' : p_mapt_solo 
                 #'centroid_solo_p': p_cent_solo,
                 #'edge_p'    : p_edge 



if __name__ == '__main__':
    pass