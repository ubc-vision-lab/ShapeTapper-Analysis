import os
import pickle
import numpy as np
import scipy.io as sio
from ShapeObj import Shape
from map_shape_pairs import Transform, rotate_points


class ShapeIO :
    
    # Called by with list of shapes to act upon, and IO paths
    def __init__(self, in_path, out_path, shape_list) :
        # Load paths, store image path for convenience
        self.in_path  = in_path
        self.out_path = out_path
        self.img_path = os.path.join(self.in_path, "Shapes")
        
        # Test shapeList to determine if it is a list of shape pairs, of shapes, or invalid format
        self.shape_list = shape_list
        # If is a list of lists of strings of length 2, then shape_list represents shape pairs
        if all(isinstance(elem, list) and len(elem)==2 and self.__checkListType(elem) for elem in shape_list) :
            self.itertype = "Pairs"
        elif self.__checkListType(shape_list) :
            self.itertype = "Singles"
        else :
            print "Error: shape list invalid. Must be list of shape names or pairs of shape names"
            self.itertype = None
            self.shape_list = None
        
        # Load shapes in shape_list
        self.__loadShapeArray()


    # Checks whether str_list is a list of only strings
    def __checkListType(self, str_list) :
        return all(isinstance(elem, basestring) for elem in str_list)


    # Iterates through shape list and stores shape objects to be acted upon by run()
    def __loadShapeArray(self) :
        if self.shape_list is None : return

        self.shapes = []
        for s in self.shape_list :
            # If single shapes, store shape objects in a list
            if self.itertype == "Singles" :
                    shape = self.__loadShape(s)
                    self.shapes.append(shape)

            # Store shape pairs in a list of lists, as in the input list
            elif self.itertype == "Pairs" :
                shape_from = self.__loadShape(s[0])
                shape_to   = self.__loadShape(s[1])
                self.__loadTrans(shape_to) # Load transform dictionary file for shape_to
                self.shapes.append([shape_from,shape_to])


    # Creates a shape object instance, with error handling
    def __loadShape (self, img_name) :
        try :
            newShape = Shape(self.img_path, img_name)
        except (IOError) :
            newShape = None
        return newShape


    # Applies function to list of shapes, looping through patient and condition lists if given
    def run(self, func, patientList=None, condList=None, taskList=None):
        if self.shape_list is None : return

        if patientList is not None:
            for p in patientList :
                print "Patient :", p
                if condList is not None :
                    for c in condList :
                        print "Condition :", c
                        if taskList is not None :
                            for t in taskList :
                                print "Task :", t
                                self.__applyToShapes(func, p, c, t)
                        else :
                            self.__applyToShapes(func, p, c)
                else :
                    self.__applyToShapes(func, p)
        else :
            self.__applyToShapes(func)


    # Applies function to list of shapes, loading observed and generated uniform points as needed
    def __applyToShapes(self, func, patient=None, cond=None, task=None) :
        if self.shape_list is None : return

        for shape in self.shapes :

            if self.itertype == "Singles" :
                if shape is None : continue
                shape.pair_mapping = None
                # Load observed touchpoint data
                shape.observed = self.__loadObs(shape, patient, task)
                # Load uniform generated data (returns None if empty)
                shape.uniform  = self.__loadUnif(shape, patient, cond)

            elif self.itertype == "Pairs" :
                shape_from = shape[0]
                shape      = shape[1]
                if shape_from is None or shape is None : continue
                # Store name of shape_from to keep track of shape pair
                shape.pair_mapping = shape_from.name
                # Load transformed observed touchpoint data from shape_from to shape
                shape.observed = self.__loadShapeObsPair(shape_from, shape, patient, task)
                # Load uniform generated data (returns None if empty)
                shape.uniform  = self.__loadUnif(shape, patient, cond)

            # Apply function to shape; 
            # patient and cond default to None in functions which do not use them
            if patient is not None:
                if cond is not None :
                    if task is not None :
                        func(shape, self.out_path, patient, cond, task)
                    else :
                        func(shape, self.out_path, patient, cond)
                else :
                    func(shape, self.out_path, patient)
            else :
                func(shape, self.out_path)


    # Load observed touch point MAT file to shape
    def __loadObs(self, shape, patient, task=None) :
        if patient is None : return None
        
        # Generate observed MAT filename
        obs_path  = os.path.join(self.in_path, patient, "observed_touchpoints")

        if task is None :
            obs_fname = "_".join((shape.name, "Patient", patient, "aggregated_observations")) + ".mat"
        else :
            obs_fname = "_".join((shape.name, "Patient", patient, "aggregated_observations", task)) + ".mat"

        # If loadMat is successfull, then add observed data to shape object
        obs_mat = self.__loadMat(obs_path, obs_fname)
        if obs_mat is None :
            if task is None :
                obs_fname = "_".join((shape.name, "Patient", patient, "observed_touchpoints")) + ".mat"
            else :
                obs_fname = "_".join((shape.name, "Patient", patient, "observed_touchpoints", task)) + ".mat"
            obs_mat = self.__loadMat(obs_path, obs_fname)
            
        if obs_mat is not None :
            observed = np.ascontiguousarray(obs_mat['img_dataset'], dtype=np.float32)
            if observed.shape[0] == 0 : return None
            observed[:,1] = shape.dims[0]-observed[:,1] # opencv coordinates use inverted y-axis
        else : observed = None
        return observed


    # Load generated uniform point MAT file to shape
    def __loadUnif(self, shape, patient, cond) :
        if patient is None or cond is None : return None
        
        # Generate uniform MAT filename
        unif_path  = os.path.join(self.in_path, patient, "uniform_points", cond)
        unif_fname = "_".join((shape.name, "Patient", patient, "uniform_points", cond)) + ".mat"

        # If loadMat is successfull, then add uniform data to shape object
        unif_mat = self.__loadMat(unif_path, unif_fname)
        if unif_mat is not None :
            uniform = np.ascontiguousarray(unif_mat['unif_datasets'], dtype=np.float32)
        else : uniform = None
        return uniform


    # Load MAT file, given its path and filename, handles read and IO errors
    def __loadMat(self, path, fname) :
        full_path = os.path.join(path, fname)
        try:
            mat = sio.loadmat(full_path)
        except (TypeError, IOError) :
            # print "Error loading: {0} -- skipping...".format(full_path)
            mat = None
        return mat


    # Loads observed points from "shape_from", fits them to "shape_to" and returns the transformed points
    def __loadShapeObsPair(self, shape_from, shape_to, patient=None, task=None) :
        if patient is None : return None
        
        # Load observed touch points from shape_from
        obs_from = self.__loadObs(shape_from, patient, task)
        if obs_from is None:
            print "Error: {0} has no observed points for patient {1}".format(shape_from.name, patient)
            return
        
        # Check shape_to transform dict to see if the matching has not been performed for shape_from already
        if shape_from.name not in shape_to.fitted_transforms :
            shape_to.fitMedialAxisFrom(shape_from) # calculate transform
            self.__saveTrans(shape_to) # save updated transform dictionary for future reference

        # Load transform parameters for shape_from -> shape_to
        tf = shape_to.fitted_transforms[shape_from.name] 
        return rotate_points(obs_from, tf, shape_from.dims) # return transformed observed points


    # Load transform dictionary from .tf file (custom binary format)
    def __loadTrans(self, shape) :
        # Generate transform dictionary file name
        tf_path  = os.path.join(self.img_path, "shape_analysis")
        tf_fname = "_".join((shape.name, "shape_pairs")) + ".tf"
        tf_full_path = os.path.join(tf_path, tf_fname)
        
        # Load transform dictionary and save in shape
        try :
            tf = pickle.load( open( tf_full_path, "rb" ) )
        except (IOError, TypeError) :
            return
        print "Loaded {0}, {1}".format(shape.name, tf)
        shape.fitted_transforms = tf


    # Save transform dictionary as .tf file (custom binary format)
    def __saveTrans(self, shape) :
        # Save transform dictionary to .tf file
        tf_path  = os.path.join(self.img_path, "shape_analysis")
        tf_fname = "_".join((shape.name, "shape_pairs")) + ".tf"
        tf_full_path = os.path.join(tf_path, tf_fname)
        pickle.dump( shape.fitted_transforms, open( tf_full_path , "wb" ) )
        
        # Save as MAT to be visible in MATLAB, for good measure
        tf_fname_mat = "_".join((shape.name, "shape_pairs")) + ".mat"
        tf_full_mat_path = os.path.join(tf_path, tf_fname_mat)
        sio.savemat(tf_full_mat_path, shape.fitted_transforms)