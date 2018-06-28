#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""Main batch processing script to run multiple steps of medial axis pipeline

@author: Jamie Dunkle
"""

from ShapeIO import ShapeIO
from generate_uniform_points import generateUniformData
from spatial_analysis  import spatialAnalysis
from generate_figures  import plotHeatMap, plotRefObjects, plotMappedShapes
from gen_cdf_figures   import plotCDFFig
from bayesian_analysis import bayesianAnalysis
 
################## Globals - CHANGE THESE #######################################
# img_names = ["blake_01","blake_03","blake_04","blake_06","blake_07",
#              "blake_08","blake_09","blake_10","blake_11","blake_12"]

### DF pairs
# img_names  = [ ["blake_04","blake_07"],
#                ["blake_07","blake_04"],
#                ["blake_11","blake_07"],
#                ["blake_07","blake_11"],
#                ["blake_10","blake_07"],
#                ["blake_07","blake_10"],
#                ["blake_04","blake_10"],
#                ["blake_10","blake_04"] ]

### MC pairs
# img_names  = [ ["blake_03","blake_07"],
#                ["blake_07","blake_03"],
#                ["blake_04","blake_07"],
#                ["blake_07","blake_04"],
#                ["blake_10","blake_07"],
#                ["blake_07","blake_10"],
#                ["blake_04","blake_10"],
#                ["blake_10","blake_04"]]

### Control pairs
# img_names  = [ ["blake_11","blake_07"],
#                ["blake_07","blake_11"],
#                ["blake_10","blake_07"],
#                ["blake_07","blake_10"],
#                ["blake_04","blake_10"],
#                ["blake_10","blake_04"]]

img_names  = [ ["blake_04","blake_07"],
               ["blake_07","blake_04"],
               ["blake_03","blake_07"],
               ["blake_07","blake_03"]]

# conditions = ["in_shape", "bounding_circle", "patient_fitted", "normal_distribution", "touchpoint_hull"]
conditions = ["in_shape", "bounding_circle", "touchpoint_hull"]

# patients   = ["DF","MC"] #"MC1" "MC2"
patients = [ 'S01','S02','S03','S04','S06','S07',
             'S08','S09','S11','S13','S14','S17',
             'S20','S21','S22','S23','S24',
             'S25','S26','S27','S28','S29','S30','S31','S32','S33']

in_path  = "D:\\ShapeTapper-Analysis\\" # where shape info and subject touch points are held
out_path = "D:\\ShapeTapper-Analysis\\" # where data and analysis files should be output


if __name__ == '__main__':

    # Load shape list, patient list and condition list into ShapeIO iterator class
    shapes = ShapeIO(in_path, out_path, img_names)


    ### RUN ANY NUMBER OF THESE FUNCTIONS BY COMMENTING/UNCOMMENTING ########
    ### They are presented in the correct order here, do not change order ###

    # shapes.run(plotRefObjects)

    # shapes.run(plotHeatMap, patients)

    # shapes.run(plotMappedShapes, patients)

    # shapes.run(plotCDFFig, patients)

    # shapes.run(generateUniformData, patients, conditions)

    shapes.run(spatialAnalysis, patients, conditions) ## relies on generateUniformData

    # shapes.run(bayesianAnalysis, patients)

