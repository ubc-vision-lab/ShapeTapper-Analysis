# ShapeTapper-Analysis
Code to analyze subject results for ShapeTapper

How to use this analysis package:

 - Make sure that you have Python 2.7 installed on your machine, with the following packages in your environment: OpenCV2 (cv2), NumPy, SciPy, Numba (instructions for installing OpenCV here: https://pypi.org/project/opencv-python/ , use `pip install opencv-python`)

 - Create an analysis folder, make sure that `in_path` and `out_path` in `ShapeTapper_Analysis.py` point here.

 - Ensure that the shape images are all in a subdirectory of the analysis folder, named `Shapes/`. You can find the current Blake shape    data in the subfolder `"Shapes/"` in this repository.

 - Put the shape analysis data in `Shapes/shape_analysis`, ex. `Shapes/shape_analysis/blake_01_shape_analysis.mat`
   (if you don't have these files, or need to generate them for a new shape, use `find_medial_axis.py`). This is also included for Blake shapes in this repository.

 - Make a patient subdirectory in this folder, for example `DF/`. [Note, you may also use preproc/parseData_from_struct.m or preproc/parseData_from_struct_Jan062020.m]
 
 - Make the following folder in the patient subdirectory: `[my_patient]/observed_touchpoints`

 - Put the touchpoint data in `[my_patient]/observed_touchpoints`, formatted as .mat files with one field (`"img_dataset"`) containing __Nx2 array__ [N points with (x,y) as columns]. 
   `preproc/` contains MATLAB scripts which should do this for you (CalculateEventInfo.m, ProcessSTFiles.m, SlimProcessSTFiles, parseData.m)

 - Open `ShapeTapper_Analysis.py` and ensure that the shape list, patient names and in/out paths defined at the top of the file are correct for your analysis.

	Note: Shape lists of pairs will map the touchpoints from the first shape to the second shape. The second shape is what will then be analysed or plotted.
	Note: Always ensure that the patient and condition lists are lists, even with one entry. I.e., use `["MC"]`, not `"MC"`

 - Comment out the functions which you do not wish to run, leave only the ones you want.

 - Open a terminal, `cd` to the directory containing your scripts and patient folder, and run `python ShapeTapper_Analysis.py`
