#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Finds best fit of medial axis point for each shape pair

@author: Jamie Dunkle
"""

import cv2
import os, errno, sys
import numpy as np
import scipy.io as sio
import distance_lib as dist
# from timeit import default_timer as timer


class Transform :
    def __init__(self, degree, offset, scale) :
        self.degree = degree
        self.offset = offset
        self.scale  = scale


################# Function Definitions #########################################################################

def rotate_points(points, tf, dims) :
    if points.shape[0]==0 : 
        return None

    # Create rotation matrix from transform object    
    rotation = cv2.getRotationMatrix2D((0,0), tf.degree, tf.scale).astype(np.float32)

    # Format points for matrix multiplication
    points = np.ascontiguousarray(points - np.array([dims[1]/2, dims[0]/2]), dtype=np.float32) # center to (0,0)
    z = np.ones((points.shape[0],1), dtype=np.float32) # need a column of ones for affine transform (2x3)
    points = np.append(points, z, axis=1) # append ones as z-axis

    # Rotate points via matrix multiplication
    points_rot = np.matmul(rotation, points.transpose()).transpose(1,0) 
    
    return (points_rot + tf.offset) # add offset


# Takes two Reference Objects (RO), usually medial axis, and find best fit via rotation & translation
def fit_ro(img_from_dims, img_from_ro, img_to_dims, img_to_ro, usefast=1) :
    
    # Get scaling factor (by canvas length) to scale shape_from RO to the size of shape_to RO
    len_from = max(img_from_dims) 
    len_to   = max(img_to_dims)
    scale    = np.true_divide(len_to, len_from)

    # Create scaled rotation matrices for angles 0 to 360
    rots = range(360)
    rotations   = [cv2.getRotationMatrix2D((0,0), r, scale) for r in rots]
    ro_centered = img_from_ro - np.array([img_from_dims[1]/2, img_from_dims[0]/2])
    z = np.ones((img_from_ro.shape[0],1), dtype=np.float32) # need a column of ones for affine transform (2x3)
    ro_centered = np.append(ro_centered, z, axis=1)

    # Apply rotation matrices to Reference Object, to create an array of rotated ROs
    ro_rots = np.array([np.matmul(r, ro_centered.transpose()) for r in rotations]).transpose(0,2,1)
    ro_rots = ro_rots.astype(np.float32)
    
    # Create 2D array of offsets
    x,y     = np.meshgrid(np.arange(img_to_dims[1]),np.arange(img_to_dims[0]))
    offsets = np.array([x,y]).transpose(1,2,0)
    offsets = offsets.astype(np.float32)

    # Limit offsets to the interquartile (1/4:3/4) range around the canvas center
    xlen = img_to_dims[0]
    ylen = img_to_dims[1]
    offsets_c = offsets[ (xlen/4):(3*xlen/4), (ylen/4):(3*ylen/4) ]
    xlen_c = offsets_c.shape[0]
    ylen_c = offsets_c.shape[1]
    
    if usefast==1 : # TWO PASS METHOD - FASTER AND USUALLY PRODUCES THE SAME RESULT AS ONE PASS, BUT PRONE TO ERRORS #########
        # First pass - test every offset in interior of shape_to, iterate by increments of 10 degrees
        inc = 10
        ro2ro_1stpass = np.empty((xlen_c, ylen_c, 360/inc)) # Array of distance measures
        for i in range(xlen_c) :
            sys.stdout.write( "\r{0} / {1}".format(i+1, xlen_c) ) # print progress counter
            sys.stdout.flush()
            for j in range(ylen_c):
                # Flatten rotated ROs into one-dimensional vector of points to speed up distance calculation
                ro_from_rot = ro_rots[range(0,360,inc)].reshape(-1, ro_rots.shape[-1])
                # Add offset
                ro_from_rot += offsets_c[i,j]
                # Find distance from rotated RO to target RO (img_to_ro)
                dist_to_ros = dist.points2points(ro_from_rot, img_to_ro)
                # Reshape flattened vector to original matrix of rotated ROs 
                dist_to_ros = dist_to_ros.reshape((360/inc,ro_rots.shape[1]))
                # Store sum min distances from RO to RO for each rotation
                ro2ro_1stpass[i,j,:] = np.sum(dist_to_ros, axis=1) 

        # We now have the best offset and approximate angle (+/-10 degrees)
        idx = np.unravel_index(np.argmin(ro2ro_1stpass, axis=None), ro2ro_1stpass.shape) # get index of min distance

        # Second pass - iterate one degree at a time to find best angle
        ro2ro_2ndpass = np.empty(2*inc)
        for i in range(2*inc) :
            ro_from_rot = ro_rots[idx[2]*inc+(i-inc)]+offsets[idx[0:2]] # from approx angle, calc dist from -10 to +10 degrees
            ro2ro_2ndpass[i] = np.sum(dist.points2points(ro_from_rot, img_to_ro))
        
        offset_min = offsets_c[idx[0:2]]
        degree_min = range(0,360)[np.argmin(ro2ro_2ndpass) - inc + (idx[2]*inc)] # unwrap min angle to get degree rotation
        ######################################################################################################################

    else : # ONE PASS METHOD - MORE ACCURATE BUT MUCH MUCH SLOWER ############################################################
        ro2ro = np.empty((xlen_c, ylen_c, 360)) # Array of distance measures
        for i in range(xlen_c) :
            sys.stdout.write( "\r{0} / {1}".format(i+1, xlen_c) ) # print progress counter
            sys.stdout.flush()
            for j in range(ylen_c) :
                # Flatten rotated ROs into one-dimensional vector of points to speed up distance calculation
                ro_from_rot  = ro_rots.reshape(-1, ro_rots.shape[-1])
                # Add offset
                ro_from_rot += offsets_c[i,j]
                # Find distance from rotated RO to target RO (img_to_ro)
                ro2ro_dist   = dist.points2points(ro_from_rot, img_to_ro)
                # Reshape flattened vector to original matrix of rotated ROs, store sum min distances
                ro2ro_dist = ro2ro_dist.reshape(ro_rots.shape[0:2])
                # Store sum min distances from RO to RO for each rotation
                ro2ro[i,j,:] = np.sum(ro2ro_dist, axis=1) 

        idx = np.unravel_index(np.argmin(ro2ro, axis=None), ro2ro.shape) # get index of min distance
        offset_min = offsets_c[idx[0:2]] # Retrieve original offset (from interquartile offset index)
        degree_min = idx[2]
        ######################################################################################################################
    
    sys.stdout.write( "\n" ) # print new line to end progress counter
    sys.stdout.flush()    
    
    tf_out = Transform(degree_min, offset_min, scale) # output a Transform object containing parameters for rotate_points()
    return tf_out



if __name__ == '__main__':
    pass