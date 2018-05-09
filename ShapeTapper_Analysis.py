from ShapeIO import ShapeIO
from generate_uniform_points import generateUniformData
from spatial_analysis import spatialAnalysis
from generate_figures import plotHeatMap, plotMedialAxis, plotMappedShapes
from gen_cdf_figures import plotCDFFig
 
################## Globals - CHANGE THESE #######################################
# img_names  = ["solo3","solo5","solo6","solo7","solo9","solo10","solo11","solo12",
#              "blake_01","blake_03","blake_04","blake_06","blake_07",
#              "blake_08","blake_09","blake_10","blake_11","blake_12"]
img_names  = [ ["blake_04","blake_11"], 
               ["blake_06","blake_11"]]

conditions = ["bounding_circle","in_shape","touchpoint_hull","patient_fitted"]
patients   = ["DF","MC"]

in_path  = "D:\\ShapeTapper-Analysis\\"
out_path = "D:\\ShapeTapper-Analysis\\"


if __name__ == '__main__':

    shapes = ShapeIO(in_path, out_path, img_names)

    # shapes.run(plotMedialAxis)

    # shapes.run(plotHeatMap, patients)

    shapes.run(plotMappedShapes, patients)

    # shapes.run(plotCDFFig, patients)

    # shapes.run(generateUniformData, patients, conditions)

    # shapes.run(spatialAnalysis, patients, conditions)
