# Image-Segmentation For Open Loop Experiment: OTSU method ( Clustered cell )
This is basic work on Image segmentation in MATLAB using OTSU method. In order to run the code, you have download the `exp_data_2` folder and copy the directory address or the path in the code. The code should be able to run in any versions of MATLAB but to be specific, I ran this code in `MATLAB2023b` version. You must also install the `Image Processing Toolbox` from the MATLAB app store. Sample results are given in file names `pos3 - new exp 2` to `pos15 - new exp 2`. To find which position in the microfluidic device has healthy cells, check the microsoft word file `2021_03_13 2iL.doc` in `exp_data_2/ON_2ILIF`.  
Now let's dive into the workflow or the function of this code:

- The main code to run is named `MAIN_OFFLINE_run.m`
- The aim to study from the image is based on clusted cell and not tracking single cell. So, the `tracking_flag = 0;` will always be same throughout the runtime of this code.
- The channel numbers represent the microscope channel input but the number is `-1` from the microscope.
- The GFP is the studied fluorescence here. All the calculated results can be found in `outputArg1.mat` file from the MATLAB workspace.
- In this particular experiment, even though Sulphorodamine channel is present but no Sulphorodamine was used in the experiment.

Now let's explore the `OFFLINE_function.m` function:
- The parameters you have entered in the MAIN section are passed here in this function.
- To note: we are not using anti-masking here, so we always select our background as Dial-A-Wave junction or any other junction available in the experiment channel 1.
- The Phase Contrast part of this code has parameters exported from the Image Segmentation app where the Cell shapes are identified as circles, dilated and eroded accordingly.
- `cb = 3` is the caliberation phase which is usually 3hrs of the experiment.


This code closely follows the original ChipSeg code from [here](https://github.com/LM-group/ChipSeg).
