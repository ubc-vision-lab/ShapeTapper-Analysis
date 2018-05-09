import os
import cv2
import scipy.io as sio
import numpy as np
from map_shape_pairs import fit_ro

class Shape :

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
            self.img[(self.img[:,:,3]==0),0:3] = 0

    def __loadRefObjs(self) :
        mat_path = os.path.join(self.img_path,'shape_analysis')
        mat_name = self.name+'_shape_analysis.mat'
        mat_full_path = os.path.join(mat_path, mat_name)
        try:
            s_mat = sio.loadmat(mat_full_path)
        except (TypeError, IOError) :
            print "Error loading: {0} -- skipping {1}...".format(mat_full_path, self.name)
            raise IOError
        self.medial_axis = np.ascontiguousarray(s_mat['ma_points'], dtype=np.float32)   # (x,y)
        self.edge_points = np.ascontiguousarray(s_mat['edge_points'], dtype=np.float32) # (x,y)
        self.centroid    = np.ascontiguousarray(s_mat['centroid'], dtype=np.float32)    # (x,y)

    def fitMedialAxisFrom(self, shapeFrom) :
        transform = fit_ro(shapeFrom.dims, shapeFrom.medial_axis, self.dims, self.medial_axis, usefast=1)
        self.fitted_transforms.update({shapeFrom.name : transform})

    def __init__(self, img_path, img_name) :
        self.img_path = img_path
        self.name = img_name
        self.fitted_transforms = dict() # stores transform objects for shape fitting
        try : 
            self.__loadImg()
            self.__loadRefObjs()
        except :
            raise IOError