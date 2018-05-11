from ShapeIO import ShapeIO
from generate_uniform_points import generateUniformData
from spatial_analysis import spatialAnalysis
from generate_figures import plotHeatMap, plotMedialAxis, plotMappedShapes
from gen_cdf_figures import plotCDFFig
 
################## Globals - CHANGE THESE #######################################
img_names = ["blake_01","blake_03","blake_04","blake_06","blake_07",
             "blake_08","blake_09","blake_10","blake_11","blake_12"]

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

conditions = ["in_shape", "bounding_circle", "patient_fitted", "normal_distribution", "touchpoint_hull"]
# conditions = ["in_shape"]

patients   = ["DF","MC"]
# patients   = ["DF"]
# patients   = ["MC"]

in_path  = "D:\\ShapeTapper-Analysis\\"
out_path = "D:\\ShapeTapper-Analysis\\"


if __name__ == '__main__':

    shapes = ShapeIO(in_path, out_path, img_names)

    # shapes.run(plotMedialAxis)

    # shapes.run(plotHeatMap, patients)

    # shapes.run(plotMappedShapes, patients)

    # shapes.run(plotCDFFig, patients)

    # shapes.run(generateUniformData, patients, conditions)

    # shapes.run(spatialAnalysis, patients, conditions)
