#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 17:02:22 2018

Finds best fit of medial axis point for each shape pair

Returns offset, rotation and scale, and new observed file

@author: Jamie Dunkle
"""

import cv2
import os, errno
import numpy as np
import scipy.io as sio
import distance_lib as dist
from timeit import default_timer as timer


class Transform :
    def __init__(self, degree, offset, scale) :
        self.degree = degree
        self.offset = offset
        self.scale  = scale


################# Function Definitions #########################################################################
def rotate_points(points, tf, dims) :
    if len(points)==0 : return None
    rotation = cv2.getRotationMatrix2D((0,0), tf.degree, tf.scale)
    points = np.ascontiguousarray(points - np.array([dims[1]/2, dims[0]/2]), dtype=np.float32)
    z = np.ones((points.shape[0],1), dtype=np.float32) # need a column of ones for affine transform
    points = np.append(points, z, axis=1)
    points_rot = np.array(np.matmul(rotation, points.transpose())).transpose(1,0)
    return points_rot + np.flip(tf.offset,axis=0)


def fit_ro(img_from_dims, img_from_ro, img_to_dims, img_to_ro, usefast=1) :
    
    # Get scaling factor (by area) to scale shape_from to the size of shape_to
    area_from = img_from_dims[0] * img_from_dims[1] 
    area_to   = img_to_dims[0] * img_to_dims[1] 
    scale     = np.true_divide(area_to, area_from)
    # len_from = max(img_from_dims) 
    # len_to   = max(img_to_dims)
    # scale    = np.true_divide(len_to, len_from)

    # Create scaled rotation matrices for angles 0:360, apply to medial axis of shape_from
    rots = range(360)
    rotations   = [cv2.getRotationMatrix2D((0,0), r, scale) for r in rots]
    ro_centered = img_from_ro - np.array([img_from_dims[1]/2, img_from_dims[0]/2])
    z = np.ones((img_from_ro.shape[0],1), dtype=np.float32) # need a column of ones for affine transform
    ro_centered = np.append(ro_centered, z, axis=1)
    ro_rots     = np.array([np.matmul(r, ro_centered.transpose()) for r in rotations]).transpose(0,2,1)
    
    # Create 2D array of offsets, take only those in the interior : (1/4:3/4) range
    x,y     = np.meshgrid(np.arange(img_to_dims[1]),np.arange(img_to_dims[0]))
    offsets = np.array([x,y]).transpose(1,2,0)
    offsets = offsets.astype(np.float32)
    ro_rots = ro_rots.astype(np.float32)
    
    if usefast==1 :
        # TWO PASS METHOD - FASTER AND USUALLY SAME AS ONE PASS, BUT PRONE TO SMALL ERRORS (1-2 DEGREES/PIXELS)
        # First pass - test every offset in interior of shape_to, iterate by increments of 10 degrees
        inc = 10
        ro2ro_1stpass = np.empty((img_to_dims[0], img_to_dims[1], 360/inc))
        for i in range(img_to_dims[0]) :
            print i+1, '/', img_to_dims[0]
            for j in range(img_to_dims[1]):
                ro_from_rot = ro_rots[range(0,360,inc)].reshape(-1, ro_rots.shape[-1]) + offsets[i,j]
                dist_to_ros = dist.points2points(ro_from_rot, img_to_ro).reshape((360/inc,ro_rots.shape[1]))
                ro2ro_1stpass[i,j,:] = np.sum(dist_to_ros, axis=1) # get min distances from RO to RO for each rotation
                    
        # Find best offset and angle (~10 degrees), iterate one degree at a time to find best angle
        idx = np.unravel_index(np.argmin(ro2ro_1stpass, axis=None), ro2ro_1stpass.shape)
        ro2ro_2ndpass = np.empty(2*inc)
        for i in range(2*inc) :
            ro_from_rot = ro_rots[idx[2]*inc+(i-inc)]+offsets[idx[0:2]]
            ro2ro_2ndpass[i] = np.sum(dist.points2points(ro_from_rot, img_to_ro))
        
        offset_min = idx[0:2]
        degree_min = range(0,360)[np.argmin(ro2ro_2ndpass) - inc + (idx[2]*inc)] 
        #######################################################################################################
    else :
        # ONE PASS METHOD - POTENTIALLY MORE ACCURATE BUT MUCH MUCH SLOWER ####################################
#        ro2ro = np.empty((img_to_dims[1], img_to_dims[0], 360))
        ro2ro = np.empty((img_to_dims[1], img_to_dims[0], 360))
        for i in range(img_to_dims[1]) :
            print i+1, '/', img_to_dims[1]
            for j in range(img_to_dims[0]) :
                ro_from_rot  = ro_rots.reshape(-1, ro_rots.shape[-1]) + offsets[i,j] # flatten all rotated ROs
                ro2ro[i,j,:] = np.sum(dist.points2points(ro_from_rot, img_to_ro).reshape(ro_rots.shape[0:2]), axis=1)
                    
        # Find best offset and angle (~10 degrees), iterate one degree at a time to find best angle
        idx = np.unravel_index(np.argmin(ro2ro, axis=None), ro2ro.shape)
        offset_min = idx[0:2]
        degree_min = idx[2]
        ######################################################################################################
    
    tf_out = Transform(degree_min, offset_min, scale)
    return tf_out