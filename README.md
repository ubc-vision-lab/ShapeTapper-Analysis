# ShapeTapper-Analysis
Code to analyze subject results for ShapeTapper

How to use this analysis package:

 - Ensure that the shape images are all in a subdirectory of this folder named `Shape/`

 - Make a patient subdirectory in this folder, for example `DF/`. 
 - Make the following folders in the patient subdirectory: `[my_patient]/aggregated_observations`, `[my_patient]/shape_analysis`

 - Put the shape analysis data in `[my_patient]/shape_analysis`, ex. `DF/shape_analysis/blake_01_shape_analysis.mat`
   (if you don't have these files, or need to generate them for a new shape, use `find_medial_axis.py`

 - Put the touchpoint data in `[my_patient]/aggregated_observations`, formatted as .mat files with one field (`"img_dataset"`) containing __Nx2 array__ [N points with (x,y) as columns]. `preproc/Transformerz.m` should do this for you.

 - Open the scripts and ensure that the patient variable defined at the top of the file matches the current patient.

 - Run the following scripts in this order:
   1. `get_dists_c_build.py` - this will build an optimized C library necessary to run the spatial analysis in a reasonable timeframe
   2. `generate_uniform_points.py` - this will generate uniform data sets in the the `generated_uniform_data/` subdirectory
   3. `spat_analysis.py` - this will perform distance and spatial analysis on the touch data and save them in `distance_analysis/` (note: this script currently takes several hours to run)
   4. `StatAnalysis.m`  - creates a statistical summary of observed touchpoints based on the uniform data sets, saved as Excel spreadsheets in the `distance_analysis/` subdirectory
   5. `heat_maps.py` - (can be run at any point) 

 Notes:
  - Most scripts have customizable parameters at the top, along with the shape names and patient name.
  - These scripts will automatically generate the subdirectories required for subsequent analysis scripts. Since these are currently hard-coded in the scripts, it is important to not change the directory structure.
