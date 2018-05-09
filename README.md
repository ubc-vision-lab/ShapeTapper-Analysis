# ShapeTapper-Analysis
Code to analyze subject results for ShapeTapper

How to use this analysis package:

 - Ensure that the shape images are all in a subdirectory of this folder named `Shape/`

 - Make a patient subdirectory in this folder, for example `DF/`. 
 - Make the following folder in the patient subdirectory: `[my_patient]/observed_touchpoints`

 - Put the shape analysis data in `Shapes/shape_analysis`, ex. `Shapes/shape_analysis/blake_01_shape_analysis.mat`
   (if you don't have these files, or need to generate them for a new shape, use `find_medial_axis.py`)

 - Put the touchpoint data in `[my_patient]/observed_touchpoints`, formatted as .mat files with one field (`"img_dataset"`) containing __Nx2 array__ [N points with (x,y) as columns]. 
   `preproc/` contains MATLAB scripts which should do this for you (CalculateEventInfo.m, ProcessSTFiles.m, SlimProcessSTFiles, parseData.m)

 - Open `ShapeTapper_Analysis.py` and ensure that the shape list, patient names and in/out paths defined at the top of the file are correct for your analysis.

	Note: Shape lists of pairs will map the touchpoints from the first shape to the second shape. The second shape is what will then be analysed or plotted.
	Note: Always ensure that the patient and condition lists are lists, even with one entry. I.e., use `["MC"]`, not `"MC"`

 - Comment out the functions which you do not wish to run, leave only the ones you want.

 - Open a terminal, `cd` to the directory containing your scripts and patient folder, and run `python ShapeTapper_Analysis.py`