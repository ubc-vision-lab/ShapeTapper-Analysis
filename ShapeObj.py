import cv2
import os
import scipy.io as sio
import numpy as np
from map_shape_pairs import fit_ro

# the Shape class stores a shape's PNG image, its name, its reference objects,
# and provides an interface to map to another shape using the map_shape_pairs library
class Shape :

    # load shape's PNG, store dimensions
    def __loadImg(self) :
        # Define full image path
        img_fname = self.name + ".png"
        img_full_path = os.path.join(self.img_path, img_fname)
        # Load image file, store in Shape object
        self.img = cv2.imread(img_full_path, cv2.IMREAD_UNCHANGED)
        if self.img is None : 
            print "Error loading: {0} -- skipping {1}...".format(img_full_path, self.name)
            raise IOError
        else :
            self.dims = self.img.shape
            self.img[(self.img[:,:,3]==0),0:3] = 0 # convert transparent pixels to black

    # load shape's Reference Object (medial axis, edge, centroid) points as a NumPy array
    def __loadRefObjs(self) :
        mat_path = os.path.join(self.img_path,'shape_analysis')
        mat_name = "_".join((self.name, "shape_analysis")) + ".mat"
        mat_full_path = os.path.join(mat_path, mat_name)
        try:
            s_mat = sio.loadmat(mat_full_path)
        except (TypeError, IOError) :
            print "Error loading: {0} -- skipping {1}...".format(mat_full_path, self.name)
            raise IOError
        self.medial_axis = np.ascontiguousarray(s_mat['ma_points'], dtype=np.float32)   # (x,y)
        self.edge_points = np.ascontiguousarray(s_mat['edge_points'], dtype=np.float32) # (x,y)
        self.centroid    = np.ascontiguousarray(s_mat['centroid'], dtype=np.float32)    # (x,y)

    # calculate the best fit of medial axis to another shape object
    # store as a transform object (rotation, offset, scale) in member dict
    def fitMedialAxisFrom(self, shape_from) :
        print "Fitting medial axis from {0} to {1}...".format(shape_from.name, self.name)
        transform = fit_ro(shape_from.dims, shape_from.medial_axis, self.dims, self.medial_axis, usefast=0)
        self.fitted_transforms.update({shape_from.name : transform}) # add to dict


    def __init__(self, img_path, img_name) :
        self.img_path = img_path
        self.name = img_name
        self.fitted_transforms = dict() # stores transform objects for shape fitting
        try : 
            self.__loadImg()
            self.__loadRefObjs()
        except :
            raise IOError