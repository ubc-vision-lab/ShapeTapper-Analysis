#!/usr/bin/python
 
import os
import cv2
import numpy as np
import scipy.io as sio
import distances as dist
import skeletonize as skel
#from timeit import default_timer as timer
    

if __name__ == '__main__':

    # Pixel offset to assist contour detection
    offset = 10

    # Contour Approximation threshold (maximum error), in pixels
    c_thresh = 2
    
    # Branch pruning threshold: >1 prunes longer branches, <1 prunes shorter
    s_thresh = 1.4

    thin2 = False # Determines whether to perform a second pass of pruning
    s_thresh2 = 0.2 # Second pass pruning threshold
    
    # Get file directories
    in_path = "Blake/" # "Simple/"
    img_names = ["blake_01","blake_04","blake_06","blake_07","blake_08","blake_09","blake_10","blake_11","blake_12"] #,"solo10","solo11","solo3","solo6","solo7","solo9"]
    out_path = './shape_analysis/'    
    
    for img_name in img_names:
        print "Starting", img_name
        
        # Read in the image.
        img_path = in_path + img_name
        img = cv2.imread(img_path,cv2.IMREAD_UNCHANGED)
        img[(img[:,:,3]==0),0:3] = 0 # Convert alpha transparency to black
        img_orig = img.copy(); # Keep a copy around
        
        # Image must be offset to find contours
        img = cv2.copyMakeBorder(img, top=offset, bottom=offset, left=offset, right=offset, 
                                 borderType=cv2.BORDER_CONSTANT, value=0);
        
        # Binarize image
        img_bin = cv2.cvtColor(img,cv2.COLOR_BGR2GRAY)
        img_bin[(img_bin!=0)] = 1 # convert to logical array (white = 1)
        
        # Compute image centroid
        white_count = np.sum(img_bin)
        white_idx = np.where(img_bin==1) # find all indices of white pixels
        centroid_y = int(np.sum(white_idx[0])/white_count)
        centroid_x = int(np.sum(white_idx[1])/white_count)
        centroid = [centroid_y,centroid_x,0]
        
        # Get piecewise linear approximation of shape edge which is threshold-pixels close
        # Allows us to perform the Integer Medial Axis transform with few unwanted branches
        contour_lines, edge_points = skel.findApproxContour(img_bin, c_thresh)
        
        # Get list of all white points (shape interior)
        white_points = np.transpose(np.array(white_idx))

        # Calculate the feature transform for interior points
        ft_dists = dist.pts2lines(white_points,contour_lines) # distance from each white point to contour
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
            ima_points, img_ima_bin = skel.pruneSkeleton(ima_points, img_ima_bin, ft_dists, s_thresh2)

        # Get distance information for centroid and medial axis points (for statistical analysis of touch data)
        dist_edges = np.min(dist.points2points(white_points,edge_points), axis=1)
        dist_medaxis = np.min(dist.points2points(white_points,ima_points),axis=1)
        dist_cent = dist.points2point(white_points,centroid)

        # Format data to x-y coordinates (OpenCV uses y-x), remove offset
        ima_points = np.flip(ima_points, axis=1) - [offset,offset]
        edge_points = np.flip(edge_points, axis=1) - [offset,offset]
        centroid = [centroid[1]-offset, centroid[0]-offset]
        
        # Draw image with medial axis points
        for p in ima_points :
            cv2.circle( img_orig, tuple(p), 1, (255,0,0,255) ) 
        for e in edge_points :
            cv2.circle( img_orig, tuple(e), 1, (0,0,255,255) ) 
        cv2.circle( img_orig, tuple(centroid), 4, (0,255,0,255) )    
        
#       # DEBUG - show image
        cv2.namedWindow('image', cv2.WINDOW_NORMAL)
        cv2.imshow('image', img_orig)
        cv2.waitKey(0)
        cv2.destroyAllWindows()

        # Save a copy of the image
        cv2.imwrite(out_path + img_name[:len(img_name)-4] + '_ma' + img_name[len(img_name)-4:],img_orig)
        
        #Save Medial Axis, Edge Points and Centroid to MAT file, along with distance information
        sio.savemat(out_path + img_name + '_shape_analysis.mat', 
                    {'ma_points':ima_points, 'edge_points':edge_points, 'centroid':centroid, 
                     'contour_lines':contour_lines, 'feature_transform':ft_dists, 'white_points':white_points, 
                     'dist_edges':dist_edges, 'dist_medaxis':dist_medaxis, 'dist_centroid':dist_cent});


