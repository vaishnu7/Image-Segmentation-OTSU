close;
clc;


mkdir mat_figures_and_data_folder %create a new folder for images

%% Standard MAIN: TO BE COMPLETED TO ANALYSE OFFLINE EXPERIMENTS
tracking_flag = 0; % 1 for tracking experiments, otherwise 0
% if (tracking_flag ==1) 
%     % in case of single cell tracking, please check the following comments
%     Parameters.fluoeval.cell_of_interest = 1;
%     % we require to identify a label for a 
%     % single cell, in order to compute a fluorescence
%     % vector associated. This lable is generally not available
%     % before a first run of the code. The authors then
%     % suggest to have a first run without the computation of
%     % fluorescence, just to detect labels associated to single cells.
% end
%%
% ----->TO MODIFY ACCORDING THE DATA OF YOUR EXPERIMENTS

n_Channels = 4; % INSERT THE NUMBER OF CHANNELS FOR THE EXPERIMENT
Blue_dye_flag = 0; % 1 if a BLUE DYE Channel is available, 0 otherwise
IRFP_flag = 1; % 1 if a NUCLEO TAG Channel is available, 0 otherwise
crop_flag = 1; % 1 if a crop is needed to analyze images, 0 otherwise
crop_drawn = 1; % 1 if a crop is manually drawn, 0 if it's downloaded from a mat file

%%
% ----> TO MODIFY ACCORDING THE TYPE OF ANALYSIS, FOR THE USE OF MASK AND ANTIMASK 

mask_channel = 1; % Insert the channel number for mask

sulpho_channel = 3; % Insert the channel number for the purple / sulpho dye
nucleo_tag_channel=4; % Insert the channel number for nucleo tag
blue_dye_channel=0; % Insert the channel number for blue dye
fluorescence_channel = 2; % Insert the channel number that will be used
% to compute fluorescence CHANNELS FLUO (e.g. GFP)

n_timeframes=43; % Insert the number of timeframes to analyse 
% (number of hour for the experiment)

position_of_interest=11; % Insert the cell position / chamber 
% of interest for the experiment 

position_of_background = 1; % Insert the position of the background - 
% usually the chamber where the junction (e.g. DAW) can be seen is 1
% (if the background is computed as the ANTIMASK insert 0)


%% Modify this part for video name and file location to access
video_name = '2i'; % Insert the file name for the VideoMaker

% Insert the name of the folder containing the images 
dir_mat_folder ='C:\Users\gz24763\OneDrive - University of Bristol\Documents\exp_data_2\ON_2ILIF\RAWDATA\2021_03_13_2iL';
% change the folder name accordingly

% to be analysed plus a slash ( \\ or /) according to the
% path method of your operative 
% system (i.e. \Mark_and_Find_001\ for Windows, Mark_and_Find_001/ for MacOS)

dir_img = '\Mark_and_Find 001\'; 
%% Parameters association (function handle - data type)
Parameters.n_timeframes = n_timeframes;
Parameters.position_of_interest = position_of_interest;
Parameters.crop_flag = crop_flag;
Parameters.crop_drawn = crop_drawn;
Parameters.mat_figures_and_data_folder = 'mat_figures_and_data_folder';
Parameters.tracking_flag = tracking_flag;
Parameters.sulpho=sulpho_channel;

Parameters.fluoeval.n_Channels = n_Channels;
Parameters.fluoeval.position_of_background = position_of_background;
Parameters.fluoeval.Blue_dye_flag = Blue_dye_flag;
Parameters.fluoeval.IRFP_flag = IRFP_flag;
Parameters.fluoeval.mask_channel = mask_channel;
Parameters.fluoeval.nucleo_tag_channel = nucleo_tag_channel;
Parameters.fluoeval.blue_dye_channel = blue_dye_channel;
Parameters.fluoeval.fluorescence_channel = fluorescence_channel;
Parameters.fluoeval.cell_of_interest = 1;

Parameters.Video.video_name = video_name;
Parameters.Video.directory_mat_folder = dir_mat_folder;
Parameters.Video.dir_img = dir_img;

% --->INSERT THE NAME OF THE EXPERIMENT FOR THE MAT FILE
save('Parameters.mat', 'Parameters');


%% OFFLINE FUNCTION CALL - NOT TESTED ONLINE / LIVE EXPERIMENT
% experiment to do first, collect time lapse image, and then do the analysis

[outputArg1,Parameters] = OFFLINE_function(Parameters);

% --->INSERT THE NAME OF THE EXPERIMENT FOR THE MAT FILE
save('Results.mat','outputArg1', 'Parameters');
% outputArg1 = stores all information about calculated fluorescence for
