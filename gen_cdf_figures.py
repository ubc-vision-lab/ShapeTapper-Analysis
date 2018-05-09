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
import distance_lib as dist
from spatial_analysis import get_cdf
from ShapeIO import ShapeIO
from timeit import default_timer as timer


################## Globals - CHANGE THESE TO RUN ON SPECIFIC SUBJECTS AND SHAPE SETS
img_names = ["blake_01","blake_03","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"]
patients = ["MC"]#,"MC","MC2"]  
img_path = "./Shapes/"         # path containing shape images
ddrive_path_prefix = "D:/ShapeTapper-Analysis/"

NUM_CDF_STEPS = 1000


################# Function Definitions #########################################################################
def drawCDF(img, img_bin, ro_pts, edge_pts, uniform, observed, regions, r_idx) :
    
    # Find region around RO corresponding to r_idx, fill with color and mark edges
    ft_dists = dist.points2points(uniform, ro_pts)
    in_region = uniform[ft_dists < np.floor(regions[r_idx])].astype(np.int32)
    in_region = np.flip(in_region,axis=1)
    region_bin = np.zeros((img.shape[0:2]), dtype=np.uint8)
    region_bin[tuple(zip(*in_region))] = 255
    region_edges = np.squeeze(cv2.findContours(region_bin,cv2.RETR_EXTERNAL,cv2.CHAIN_APPROX_NONE)[1][0])
#    cv2.fillConvexPoly(img, region_edges, (255,0,255,100))
    ### NEED TO DRAW LINES NOT CIRCLES
#    for hp in region_edges :
#        cv2.circle( img, (int(hp[0]),int(hp[1])), 1, (255,0,255,200), thickness=-1, lineType=cv2.LINE_AA )
    
    for p in in_region :
        cv2.circle( img, (int(p[1]),int(p[0])), 1, (255,0,255,100), thickness=-1, lineType=cv2.LINE_AA )
    
    # Seperate observed touch points which are in vs out of the region    
    obs_in_region = []
    obs_outof_region = []
    for i, tp in enumerate(observed) :
        if cv2.pointPolygonTest(region_edges, tuple(tp), measureDist=False) == 1 :
            obs_in_region.append(tp)
        else :
            obs_outof_region.append(tp)
    obs_in_region = np.array(obs_in_region, dtype=np.float32)
    obs_outof_region = np.array(obs_outof_region, dtype=np.float32)
    
    num_purple = img[np.all(img[:,:,0:3]==[255,0,255], axis=2)].shape[0]
    num_white  = img[np.all(img[:,:,0:3]==[255,255,255], axis=2)].shape[0]
    
    print np.true_divide(in_region.shape[0], uniform.shape[0]), "(points)", "<-->", np.true_divide(num_purple, num_purple+num_white), "(pixels)"
    
    # Draw each set of touch points with different colors
    for op in obs_in_region :
        cv2.circle( img, (int(op[0]),int(op[1])), 2, (0,0,255,255), thickness=-1, lineType=cv2.LINE_AA )  
    for op in obs_outof_region :
        cv2.circle( img, (int(op[0]),int(op[1])), 2, (255,0,0,155), thickness=-1, lineType=cv2.LINE_AA )

    for mp in ro_pts :
        cv2.circle( img, (int(mp[0]),int(mp[1])), 1, (255,0,0,255), thickness=-1, lineType=cv2.LINE_AA )
    
    for ep in edge_pts :
        cv2.circle( img, (int(ep[0]),int(ep[1])), 1, (0,0,0,255), thickness=-1, lineType=cv2.LINE_AA )
    
    return img


def plotCDFFig(shape, out_path, patient, cond=None):    
    ################### Load data ###################

    print "Plotting {0}".format(shape.name)

    if patient is None :
        print "Error in plotHeatMap() : no patient name specified."
        return

    if shape.observed is None :
        print "Error in plotHeatMap() : no observed touch points found for {0}. Please check data files.".format(shape.name)
        return

    out_path = os.path.join(out_path, patient, "figures", "cdf_figs")
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    img = shape.img.copy()

    observed = shape.observed
    medial_axis = shape.medial_axis.astype(np.uint32)
    edge_points = shape.edge_points.astype(np.uint32)

    # Trim out of shape touchpoints for presentation
    observed_inshape = [] 
    for i in range(observed.shape[0]) :
        if cv2.pointPolygonTest(edge_points, tuple(observed[i]), measureDist=False) == 1 :
            observed_inshape.append(observed[i])
    observed = np.array(observed_inshape)
    
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
    
#    # Use a slightly enlarged regions vector for centroid to catch all points
#    max_r_ct = int(np.ceil(bd_circ_radius*np.sqrt(2))) 
#    regions_ct = np.linspace(0.0, max_r_ct, num=NUM_CDF_STEPS, endpoint=True)
#    regions_ct = np.ascontiguousarray(regions_ct, dtype=np.float32)
    
    cdfs = get_cdf(medial_axis, regions_bd, shape_points, observed)
    
    onethird_idx  = np.where(cdfs[0]>np.true_divide(1,3))[0][0]
    twothirds_idx = np.where(cdfs[0]>np.true_divide(2,3))[0][0]
    
    # Calculate spatial CDF data for each Reference Object (Medax, Edge, Centroid)
    cdf_ma_img_1 = drawCDF(img.copy(), img_bin.copy(), medial_axis, edge_points, shape_points, observed, regions_bd, onethird_idx)
    cdf_ma_img_2 = drawCDF(img.copy(), img_bin.copy(), medial_axis, edge_points, shape_points, observed, regions_bd, twothirds_idx)
#    cdf_edge = plot_cdf(edge_points, regions_bd, shape_points, observed)
#    cdf_cent = plot_cdf(centroid, regions_ct, shape_points, observed)
    
    if shape.pair_mapping is not None :
        out_fname_1 = "_".join((shape.pair_mapping, "to", shape.name, "Patient", patient, '_cdf_ma_onethird.png'))
        out_fname_2 = "_".join((shape.pair_mapping, "to", shape.name, "Patient", patient, '_cdf_ma_twothirds.png'))
    else :
        out_fname_1 = "_".join((shape.name, "Patient", patient, '_cdf_ma_onethird.png'))
        out_fname_2 = "_".join((shape.name, "Patient", patient, '_cdf_ma_twothirds.png'))

    cv2.imwrite( os.path.join(out_path,out_fname_1), cdf_ma_img_1 )
    cv2.imwrite( os.path.join(out_path,out_fname_2), cdf_ma_img_2 )


if __name__ == '__main__':

    img_names  = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
                  "blake_01","blake_03","blake_04","blake_06","blake_07",
                  "blake_08","blake_09","blake_10","blake_11","blake_12"]

    patients   = ["DF","MC"]

    in_path  = "D:/ShapeTapper-Analysis/"
    out_path = "D:/ShapeTapper-Analysis/"

    shapes = ShapeIO(in_path, out_path, img_names)

    shapes.run(plotCDFFig, patients)