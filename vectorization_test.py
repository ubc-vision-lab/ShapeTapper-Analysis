# -*- coding: utf-8 -*-
"""
Created on Wed Feb 28 13:56:47 2018

@author: Visionlab
"""
import unittest
import inspect
import numpy as np
import numpy.testing as npt
import vectorization as vect

class VectorizationTestCase(unittest.TestCase):
    
    def setUp(self):
        self.line_set1 = np.array( [ [[0,0],[1,1]], [[1,1],[2,2]] ] ) # no nulls
        self.line_set2 = np.array( [ [[0,0],[1,1]], [[1,1],[1,1]] ] ) # with nulls
        self.point_set1 = np.array( [[0,1],[0,2]] )
        self.point_set2 = np.array( [[0,1],[-1,2],[0,2]] )
        
        self.correct_distances_set1 = [np.sqrt(2)/2,np.sqrt(2)] # no nulls
        self.correct_distances_set2 = [np.sqrt(2)/2,np.sqrt(5)] # with nulls
    
    def test_correct_distances_no_nulls(self):
        self.logPoint()
        self.assertAlmostEqual(vect.dist_pts2lines(self.point_set1,self.line_set1).tolist(),self.correct_distances_set1,places=7)
        
    def test_correct_distances_with_nulls(self):
        self.logPoint()
        self.assertAlmostEqual(vect.dist_pts2lines(self.point_set2,self.line_set2).tolist(),self.correct_distances_set2,places=7)
 
    def logPoint(self):
        currentTest = self.id().split('.')[-1]
        callingFunction = inspect.stack()[1][3]
        print 'in %s - %s()' % (currentTest, callingFunction)

if __name__ == '__main__':
    #unittest.main()
    line_set1 = np.array( [ [[0,0],[1,1]], [[1,1],[2,2]] ] ) # no nulls
    line_set2 = np.array( [ [[0,0],[1,1]], [[1,1],[1,1]] ] ) # with nulls
    point_set1 = np.array( [[0,1],[0,2]] )
    point_set2 = np.array( [[0,1],[-1,2],[0,2]] )
    correct_distances_set1 = [np.sqrt(2)/2,np.sqrt(2)] # no nulls
    correct_distances_set2 = [np.sqrt(2)/2,np.sqrt(5),np.sqrt(2)] # with nulls
    npt.assert_almost_equal(vect.dist_pts2lines(point_set1,line_set1).tolist(),correct_distances_set1,decimal=100)
    npt.assert_almost_equal(vect.dist_pts2lines(point_set2,line_set2).tolist(),correct_distances_set2,decimal=100)