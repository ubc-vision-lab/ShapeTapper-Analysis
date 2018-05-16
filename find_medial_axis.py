#!/usr/bin/python
 
import os, errno
import cv2
import numpy as np
import scipy.io as sio
import distance_lib as dist
import skeletonize_lib as skel
#from timeit import default_timer as timer

offset = 10

def findMedialAxis(shape, out_path, patient=None, cond=None, c_thresh=2, s_thresh=1.2, s_thresh2=0.2, thin2=False) :

    print "Finding medial axis for {0}".format(shape.name)

    out_path = os.path.join(out_path, 'shape_analysis')
    try:
        os.makedirs(out_path)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    img_orig = shape.img.copy()

    # Image must be offset to find contours
    img = cv2.copyMakeBorder(img_orig, top=offset, bottom=offset, left=offset, right=offset, 
                                 borderType=cv2.BORDER_CONSTANT, value=0)

    # Binarize image
    img_bin = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
    img_bin[(img_bin!=0)] = 1 # convert to logical array (white = 1)
    
    # Compute image centroid
    white_count = np.sum(img_bin)
    white_idx   = np.where(img_bin==1) # find all indices of white pixels
    centroid_y = int(np.sum(white_idx[0]) / white_count)
    centroid_x = int(np.sum(white_idx[1]) / white_count)
    centroid = [centroid_y, centroid_x, 0]
    
    # Get piecewise linear approximation of shape edge which is threshold-pixels close
    # Allows us to perform the Integer Medial Axis transform with few unwanted branches
    contour_lines, edge_points = skel.findApproxContour(img_bin, c_thresh)
    
    # Get list of all white points (shape interior)
    white_points = np.transpose(np.array(white_idx))

    # Calculate the feature transform for interior points
    ft_dists = dist.pts2lines(white_points, contour_lines) # distance from each white point to contour
    ft_lines = np.argmin(ft_dists,axis=1) # contour line closest to each white point
    ft_dists = np.min(ft_dists, axis=1)

    # Find the Integer Medial Axis
    ima_points = skel.imaTransform(white_points, edge_points, ft_lines)

    # Create a binary img representation of IMA points
    img_ima_bin = np.zeros_like(img_bin)
    img_ima_bin[ima_points[:,0],ima_points[:,1]] = 1
    
    # Prune branches on the IMA caused by joints in the contour line segments
    ima_points, img_ima_bin = skel.pruneSkeleton(ima_points, white_points, img_ima_bin, ft_dists, s_thresh)

    if thin2:
        ima_points, img_ima_bin = skel.pruneSkeleton(ima_points, white_points, img_ima_bin, ft_dists, s_thresh2)

    # Format data to x-y coordinates (OpenCV uses y-x), remove offset
    ima_points = np.flip(ima_points, axis=1) - [offset,offset]
    edge_points = np.flip(edge_points, axis=1) - [offset,offset]
    centroid = [ centroid[1]-offset, centroid[0]-offset ]
    
    # Draw image with medial axis points
    for p in ima_points :
        cv2.circle( img_orig, tuple(p), 1, (255,0,0,255) ) 
    for e in edge_points :
        cv2.circle( img_orig, tuple(e), 1, (0,0,255,255) ) 
    cv2.circle( img_orig, tuple(centroid), 4, (0,255,0,255) )    
    
    # DEBUG - show image
    # cv2.namedWindow('image', cv2.WINDOW_NORMAL)
    # cv2.imshow('image', img_orig)
    # cv2.waitKey(0)
    # cv2.destroyAllWindows()

    # Save after hand cleaning img_ima_bin in Spyder
#        x = np.array(np.where(img_ima_bin==1))
#        x = np.array(zip(*x))
#        ima_points = np.flip(x, axis=1) - [offset,offset]

    # Save a copy of the image
    cv2.imwrite(out_path + shape.name + '_ma.png', img_orig)
    
    #Save Medial Axis, Edge Points and Centroid to MAT file, along with distance information
    sio.savemat(out_path + shape.name + '_shape_analysis.mat', 
                {'ma_points':ima_points, 'edge_points':edge_points, 'centroid':centroid})


if __name__ == '__main__':
    pass