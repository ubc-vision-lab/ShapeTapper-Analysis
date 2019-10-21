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
img_names = ["blake_01","blake_03","blake_04","blake_06","blake_07",
             "blake_08","blake_09","blake_10","blake_11","blake_12"]
             
# img_names = ["blake_01","blake_04","blake_07",
#              "blake_10","blake_11","blake_12"]

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

# img_names  = [ ["blake_04","blake_07"],
#                ["blake_07","blake_04"],
#                ["blake_03","blake_07"],
#                ["blake_07","blake_03"]]

# conditions = ["in_shape", "bounding_circle", "patient_fitted", "normal_distribution", "touchpoint_hull"]
conditions = ["in_shape", "bounding_circle"]

tasks = None
# tasks = ["Simultaneous_2AFC","Sequential_2AFC","Oddball","N_Back"]
# tasks = ["Simultaneous_2AFC","Sequential_2AFC","Oddball","N_Back"]
# tasks = ['LR_167ms','LR_233ms','LR_500ms','LR_10000ms','OB_167ms','OB_233ms','OB_500ms','OB_10000ms']

# patients   = ["MC1", "MC2"]
# patients   = ["DF","MC"]

# patients = ['S01','S02']
# patients = ['S03','S04']
# patients = ['S05','S06']
# patients = ['S07','S08']
# patients = ['S09','S10']
# patients = ['S11','S12','S13']

# patients = ['S01','S02','S03','S04','S05','S06']
# patients = ['S07','S08','S09','S10','S11','S12']
# patients = ['S13','S14','S15','S16','S17','S18']
# patients = ['S19','S20','S21','S22','S23','S24']
# patients = ['S25','S26','S27','S28','S29','S30']
# patients = ['S31','S32','S33','S34','S35','S36']
# patients = ['S37','S38','S39','S40','S41','S42']
# patients = ['S43','S44','S45','S46','S47','S48']
# patients = ['S49','S50','S51','S52','S53','S54']

# patients = ['S01','S02','S03','S04','S05','S06',
#             'S07','S08','S09','S10','S11','S12',
#             'S13','S14','S15','S16','S17','S18',
#             'S19','S20','S21','S22','S23','S24']
# patients = ['S25','S26','S27','S28','S29','S30',
#             'S31','S32','S33','S34','S35','S36',
#             'S37','S38','S39','S40','S41','S42',
#             'S43','S44','S45','S46','S47','S48',
#             'S49','S50','S51','S52','S53','S54']

# in_path  = "D:\\ShapeTapper-Analysis\\" # where shape info and subject touch points are held
# out_path = "D:\\ShapeTapper-Analysis\\" # where data and analysis files should be output

in_path  = "E:\\ShapeTapper-Analysis\\" # where shape info and subject touch points are held
out_path = "E:\\ShapeTapper-Analysis\\" # where data and analysis files should be output


if __name__ == '__main__':

    # Load shape list, patient list and condition list into ShapeIO iterator class
    shapes = ShapeIO(in_path, out_path, img_names)


    ### RUN ANY NUMBER OF THESE FUNCTIONS BY COMMENTING/UNCOMMENTING ########
    ### They are presented in the correct order here, do not change order ###

    # shapes.run(plotRefObjects)

    # shapes.run(plotHeatMap, patients)

    # shapes.run(plotMappedShapes, patients)

    # shapes.run(plotCDFFig, patients)

    # shapes.run(generateUniformData, patients, conditions, tasks)

    # shapes.run(spatialAnalysis, patients, conditions, tasks) ## relies on generateUniformData

    # shapes.run(bayesianAnalysis, patients)

