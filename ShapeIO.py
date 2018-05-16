import os
import pickle
import numpy as np
import scipy.io as sio
from ShapeObj import Shape
from map_shape_pairs import Transform, rotate_points

class ShapeIO :

    def __loadShape (self, img_name) :
        try :
            newShape = Shape(self.img_path, img_name)
        except (IOError) :
            newShape = None
        return newShape


    def __loadMat(self, path, fname) :
        full_path = os.path.join(path, fname)
        try:
            mat = sio.loadmat(full_path)
        except (TypeError, IOError) :
            # print "Error loading: {0} -- skipping...".format(full_path)
            mat = None
        return mat


    def __loadObs(self, shape, patient) :
        if patient is None : return None
        obs_path  = os.path.join(self.in_path, patient, "observed_touchpoints")
        obs_fname = "_".join((shape.name, "Patient", patient, "observed_touchpoints")) + ".mat"
        obs_mat = self.__loadMat(obs_path, obs_fname)
        if obs_mat is not None :
            observed = np.ascontiguousarray(obs_mat['img_dataset'], dtype=np.float32)
            if observed.shape[0] == 0 : return None
            observed[:,1] = shape.dims[0]-observed[:,1] # opencv coordinates use inverted y-axis
        else : observed = None
        return observed


    def __loadUnif(self, shape, patient, cond) :
        if patient is None or cond is None : return None
        unif_path  = os.path.join(self.in_path, patient, "uniform_points", cond)
        unif_fname = "_".join((shape.name, "Patient", patient, "uniform_points", cond)) + ".mat"
        unif_mat = self.__loadMat(unif_path, unif_fname)
        if unif_mat is not None :
            uniform = np.ascontiguousarray(unif_mat['unif_datasets'], dtype=np.float32)
        else : uniform = None
        return uniform


    def __loadShapeObsPair(self, shape_from, shape_to, patient=None) :
        if patient is None : return None
        obs_from = self.__loadObs(shape_from, patient)
        if obs_from is None:
            print "Error: {0} has no observed points for patient {1}".format(shape_from.name, patient)
            return
        if shape_from.name not in shape_to.fitted_transforms :
            shape_to.fitMedialAxisFrom(shape_from)
            self.__saveTrans(shape_to)
        tf = shape_to.fitted_transforms[shape_from.name]
        return rotate_points(obs_from, tf, shape_from.dims)


    def __loadTrans(self, shape) :
        tf_path  = os.path.join(self.img_path, "shape_analysis")
        tf_fname = "_".join((shape.name, "shape_pairs")) + ".tf"
        tf_full_path = os.path.join(tf_path, tf_fname)
        try :
            tf = pickle.load( open( tf_full_path, "rb" ) )
        except (IOError, TypeError) :
            return
        print "Loaded {0}, {1}".format(shape.name, tf)
        shape.fitted_transforms = tf


    def __saveTrans(self, shape) :
        tf_path  = os.path.join(self.img_path, "shape_analysis")
        tf_fname = "_".join((shape.name, "shape_pairs")) + ".tf"
        tf_full_path = os.path.join(tf_path, tf_fname)
        pickle.dump( shape.fitted_transforms, open( tf_full_path , "wb" ) )
        tf_fname_mat = "_".join((shape.name, "shape_pairs")) + ".mat"
        tf_full_mat_path = os.path.join(tf_path, tf_fname_mat)
        sio.savemat(tf_full_mat_path, shape.fitted_transforms)


    def __checkListType(self, strlist) :
        return all(isinstance(elem, basestring) for elem in strlist)


    def __init__(self, in_path, out_path, shapeList) :
        self.in_path  = in_path
        self.out_path = out_path
        self.img_path = os.path.join(self.in_path, "Shapes")
        
        # Test shapeList to determine if it is a list of shape pairs or shapes
        self.shape_list = shapeList
        if all(len(elem)==2 and self.__checkListType(elem) and isinstance(elem, list) for elem in shapeList) :
            self.itertype = "Pairs"
        elif self.__checkListType(shapeList) :
            self.itertype = "Singles"
        else :
            print "Error: shape list invalid. Must be list of shape names or pairs of shape names"
            self.itertype = None
            self.shape_list = None
        self.__loadShapeArray()


    def __loadShapeArray(self) :
        if self.shape_list is None : return
        self.shapes = []
        for s in self.shape_list :
            if self.itertype == "Singles" :
                    shape = self.__loadShape(s)
                    self.shapes.append(shape)

            elif self.itertype == "Pairs" :
                shape_from = self.__loadShape(s[0])
                shape_to   = self.__loadShape(s[1])
                self.__loadTrans(shape_to)
                self.shapes.append([shape_from,shape_to])


    def __applyToShapes(self, func, patient=None, cond=None) :
        if self.shape_list is None : return

        for shape in self.shapes :

            if self.itertype == "Singles" :
                if shape is None : continue
                shape.pair_mapping = None
                shape.observed = self.__loadObs(shape, patient)
                shape.uniform  = self.__loadUnif(shape, patient, cond)

            elif self.itertype == "Pairs" :
                shape_from = shape[0]
                shape      = shape[1]
                if shape_from is None or shape is None : continue
                shape.pair_mapping = shape_from.name
                shape.observed = self.__loadShapeObsPair(shape_from, shape, patient)
                shape.uniform  = self.__loadUnif(shape, patient, cond)

            if patient is not None:
                if cond is not None :
                    func(shape, self.out_path, patient, cond)
                else :
                    func(shape, self.out_path, patient)
            else :
                func(shape, self.out_path)


    def run(self, func, patientList=None, condList=None):
        if self.shape_list is None : return

        if patientList is not None:
            for p in patientList :
                print "Patient :", p
                if condList is not None :
                    for c in condList :
                        print "Condition :", c
                        self.__applyToShapes(func, p, c)
                else :
                    self.__applyToShapes(func, p)
        else :
            self.__applyToShapes(func)
