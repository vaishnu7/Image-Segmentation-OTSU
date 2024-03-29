function [Output,Parameters] = OFFLINE_function(Parameters) % function created

close
clc
warning off;

%% Definition of parameters
n_CHANNELS = Parameters.fluoeval.n_Channels;
n_timeframes = Parameters.n_timeframes;
selected_BCKGfield = Parameters.fluoeval.position_of_background; %if 0, background as antimask


%% Fluorescence vectors
Output.Absolute_Fluorescence_AVG = [];
Output.Fluorescence_Channel_A_A = [];
Output.Nuclei_Normalized_Fluorescence_Channel_A_A = [];
Output.Average_Background_fluorescence = [];
Output.MASK = {};

%% CROP FUNCTION call
Parameters = CROPfunction(Parameters);

%% MAIN LOOP with mask and fluorescence computation
for i=1:n_timeframes
    % Function defining addresses associated to images
    [~, fileNames, fileNamesBCKGRND] = fileNames_addressFunction(i, Parameters);
       
    for j= 1:n_CHANNELS % updating each channel numbers to fileNamesBCKGRND{j} 
        % e.g. 5 in this experiment
        Parameters.fluoeval.CHANNEL{j} = fileNames{j};
        if (selected_BCKGfield~=0)
            Parameters.fluoeval.BCKGRND_CHANNEL{j} = fileNamesBCKGRND{j};
        end
    end
    

        % cell clusters mask
        [mask, Output.area_mask(i)] = mask_function(Parameters.fluoeval.CHANNEL{Parameters.fluoeval.mask_channel},Parameters.fluoeval.crop);
        [Output_fluo_eval] = fluorescence_evaluation(Parameters.fluoeval, mask);
        Output.MASK{i} = mask;

    % saving fluorescence values
    Output.Absolute_Fluorescence_AVG = [Output.Absolute_Fluorescence_AVG ; Output_fluo_eval.Absolute_Fluorescence_AVG];
    Output.Fluorescence_Channel_A_A = [ Output.Fluorescence_Channel_A_A ; Output_fluo_eval.Fluorescence_Channel_A_A ];
    
    % if IRFP or Nucleo-Tag is not present, then replace this with the
    % fluorescence of interest
    if (Parameters.fluoeval.IRFP_flag~=0)
        Output.Nuclei_Normalized_Fluorescence_Channel_A_A = [ Output.Nuclei_Normalized_Fluorescence_Channel_A_A ; Output_fluo_eval.Nuclei_Normalized_Fluorescence_Channel_A_A ];
    end
    Output.Average_Background_fluorescence = [ Output.Average_Background_fluorescence ; Output_fluo_eval.Average_Background_fluorescence ];
    
    
end
% VIDEO GENERATOR
Generate_Video(Parameters,Output);

% saving current workspace
curr_filename = strcat('field_',num2str(Parameters.position_of_interest));
save(strcat('./',Parameters.mat_figures_and_data_folder,'/',curr_filename));
%  strcat: horizontally concatenates the text in its input arguments e.g. A + B = AB
end

%% FUNCTION FOR THE CROP
function [Parameters] = CROPfunction(Parameters)

% Input: stack of parameters that are set up in the main file.
% Output: stack of parameters updated with the crop variable.

source = Parameters.Video.directory_mat_folder ;
imageDir_exp = Parameters.Video.dir_img ;
position_of_interest = Parameters.position_of_interest;
selected_BCKGfield = Parameters.fluoeval.position_of_background;

if (Parameters.crop_flag==1) 
    if (Parameters.crop_drawn==1) 
        % names of images must be changed according to the experimental setup and users standards
        name_position_ctrld = strcat('Position', ...
            pad(num2str(position_of_interest), 3, 'left', '0'), ...
            '--t', pad(num2str(0), 2, 'left', '0'),...
            '--C', pad(num2str(0), 2, 'left', '0'), '.tif');
        %-----------------------------------------------------------
        image = strcat(source,imageDir_exp,name_position_ctrld);
        [~,crop] = imcrop(imadjust(imread(image)));

        % imcrop: creates an interactive crop image tool
        % ------- allows us to crop the image manually using mouse 
        % imadjust: maps the intensity values in grayscale image 
        % imread: read the image from the file specified by file name

        save('cropfiles.mat','crop'); 
        Parameters.fluoeval.crop = crop;
        if (Parameters.fluoeval.position_of_background~=0)
            % name of background image must be changed 
            % according to the experimental setup and user standards
            name_BCKGRND_position = strcat('Position', ...
                pad(num2str(selected_BCKGfield), 3, 'left', '0'), ...
                '--t', pad(num2str(0), 2, 'left', '0'),...
                '--C', pad(num2str(0), 2, 'left', '0'), '.tif');
            %-------------------------------------------------------------
            imageBCKG = strcat(source,imageDir_exp,name_BCKGRND_position);
            [~,CropRectBCKG]=imcrop(imadjust(imread(imageBCKG)));
            save('cropfiles.mat','CropRectBCKG','-append');
            Parameters.fluoeval.CropRectBCKG = CropRectBCKG;
        end
    else
        % the image we want to study
        crop_var = strcat('CropRect',num2str(position_of_interest));
        temp_struct = load(Parameters.name_matfile,crop_var);
        crop = temp_struct.(crop_var);
        Parameters.fluoeval.crop = crop;
        
        if (selected_BCKGfield~=0)
            temp_struct = load(Parameters.name_matfile,'CropRectBCKG');
            CropRectBCKG = temp_struct.('CropRectBCKG');
            Parameters.fluoeval.CropRectBCKG = CropRectBCKG;
        end
    end
else
    Parameters.fluoeval.crop = 0;
end

end

%% FUNCTION TO READ CROPPED IMAGES AND GET IMAGE PATHS

function [cropped_image,fileNames, fileNamesBCKGRND] = fileNames_addressFunction(i, Parameters)
% Input: a struct of parameters, and the current time 'i'
% Output: a n*m matrix (cropped_image), where n*m depends on the crop size,
% the complete address to the position of interest (fileNames)
% and the complete address to the background position (fileNamesBCKGRND)

% Pre-alocation of used cell arrays -- {}
% change the size of cell arrays according to the size
% of your image and no. of pixel needed
fileNames = {zeros(2,2,'uint16')};
uncropped_image = {zeros(600,600,'uint16')};
cropped_image = {zeros(500,500,'uint16')};
fileNamesBCKGRND = {zeros(2,2,'uint16')};

for type_i = 1:Parameters.fluoeval.n_Channels
    % names of images must be changed according to the experimental setup
    position_contld = strcat('Position', ...
        pad(num2str(Parameters.position_of_interest), 3, 'left', '0'), ...
        '--t', pad(num2str(i), 2, 'left', '0'),...
        '--C', pad(num2str(type_i - 1), 2, 'left', '0'), '.tif'); 
    % image is identified by the file name e.g. Position011-t01-C00.tif
    %---------------------------------------------------------------
    ctrl_position = fullfile(Parameters.Video.dir_img , position_contld);
    fileNames{type_i} = fullfile(Parameters.Video.directory_mat_folder , ctrl_position);
    uncropped_image{type_i} = imread(strcat(Parameters.Video.directory_mat_folder ,ctrl_position));
    
    if (Parameters.crop_flag~=0)
        cropped_image{type_i}=imcrop(uncropped_image{type_i},Parameters.fluoeval.crop);
    else
        cropped_image{type_i}=uncropped_image{type_i};
    end
    
    if (Parameters.fluoeval.position_of_background~=0)
        % names of images must be changed according to the experimental setup
        BCKGRND_position = strcat('Position', ...
            pad(num2str(Parameters.fluoeval.position_of_background), 3, 'left', '0'), ...
            '--t', pad(num2str(i), 2, 'left', '0'),...
            '--C', pad(num2str(type_i - 1), 2, 'left', '0'), '.tif');
        %---------------------------------------------------------------
        BCKGRND_position = fullfile(Parameters.Video.dir_img , BCKGRND_position);
        fileNamesBCKGRND{type_i} = fullfile(Parameters.Video.directory_mat_folder, BCKGRND_position);
        
    end
    
end
end


%% MASK FUNCTION - USED FOR CELL CLUSTERS EXPERIMENTS

function [mask,areamask] = mask_function(Image,crop)
% Input: image of interest (Image), depending on the channel used to compute the
% mask, and the crop size n*m (crop).
% Outputs: n*m mask matrix  (mask) and associated area ( areamask).

Image_temp1=imread(Image);
if (crop~=0)
    Image_temp2=imcrop(Image_temp1,crop);
else
    Image_temp2 = Image_temp1;
end
Img = imadjust(Image_temp2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter values to change as per needed - Check MATLAB Image Segmenter App
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PHASE CONTRAST


BW1 = uint16(Img); % uint16 = 16-bit unsigned integer - takes only non-negative values
% represents values from 0 - 65,535

I = BW1;
BW2 =I> 17430;

radius=3;
decomposition=0;
se = strel('disk', radius, decomposition);
BW = imdilate(BW2, se);

BW = imfill(BW, 'holes'); 
% fills holes in the input binary image BW. 
% In this syntax, a hole is a set of background pixels that 
% cannot be reached by filling in the background from the edge of the image

radius = 8;
decomposition = 0;
se = strel('disk', radius, decomposition);
BW = imerode(BW, se);
BWfinal=BW;
  
%% as usual

if (numel(find(BWfinal))/numel(BWfinal)>=.95)
    BWfinal = ones(size(Image));
end

mask = BWfinal;

areamask=sum(sum(mask));

end

%% VIDEO GENERATOR
function [] = Generate_Video(Parameters,Output) % MATLAB default function
profile = 'MPEG-4';

video_name = Parameters.Video.video_name ;
time_lapse = VideoWriter(strcat(video_name,'_field_',num2str(Parameters.position_of_interest)),profile);
time_lapse.FrameRate = 10; % value that modify the speed
time_lapse.Quality = 100;
open(time_lapse); % open file for writing video data
fig1 = figure;


n_CHANNELS = Parameters.fluoeval.n_Channels;
n_timeframes = Parameters.n_timeframes;
selected_BCKGfield = Parameters.fluoeval.position_of_background;
% tracking_flag = Parameters.tracking_flag;
%% Added to normalize the comand imadjust for images visualiatation
% variables allocation
min_daw=inf;
max_daw=0;
min_gfp=inf;
max_gfp=0;
min_h2b=inf;
max_h2b=0;
min_sulpho=inf;
max_sulpho=0;
min_sulpho_crop=inf;
max_sulpho_crop=0;

for i=1:n_timeframes
    [cropped_image,~, fileNamesBCKGRND] = fileNames_addressFunction(i, Parameters);
    %% Check this channel whether it's mcherry in exp or make it general with variable out
     daw=imread(fileNamesBCKGRND{Parameters.sulpho}); %%%%%% check this part
     min_daw=min([min(min(daw)), min_daw]);
     max_daw=max([max(max(daw)), max_daw]);
     
     gfp=cropped_image{Parameters.fluoeval.fluorescence_channel};
     min_gfp=min([min(min(gfp)), min_gfp]);
     max_gfp=max([max(max(gfp)), max_gfp]);
     
     h2b=cropped_image{Parameters.fluoeval.nucleo_tag_channel};
     min_h2b=min([min(min(h2b)), min_h2b]);
     max_h2b=max([max(max(h2b)), max_h2b]);
     
     sulpho=cropped_image{4};
     min_sulpho=min([min(min(sulpho)), min_sulpho]);
     max_sulpho=max([max(max(sulpho)), max_sulpho]);

     %% control channel 
     sulpho_crop=imcrop(imread(Parameters.fluoeval.BCKGRND_CHANNEL{Parameters.sulpho}),Parameters.fluoeval.CropRectBCKG);
     min_sulpho_crop=min([min(min(sulpho_crop)), min_sulpho_crop]);
     max_sulpho_crop=max([max(max(sulpho_crop)), max_sulpho_crop]);

end


for i=1:n_timeframes

    [cropped_image,fileNames, fileNamesBCKGRND] = fileNames_addressFunction(i, Parameters);
    
    %% FLUO EVAL

    % 
    % for j= 1:n_CHANNELS
    %     Parameters.fluoeval.CHANNEL{j} = fileNames{j};
    %     if (selected_BCKGfield~=0)
    %         Parameters.fluoeval.BCKGRND_CHANNEL{j} = fileNamesBCKGRND{j};
    %     end
    % end
    Parameters.fluoeval.CHANNEL = fileNames(1:n_CHANNELS);
    if (selected_BCKGfield~=0)
        Parameters.fluoeval.BCKGRND_CHANNEL = fileNamesBCKGRND(1:n_CHANNELS);
    end

%% PLOTS IMAGES of interest TO MODIFY 
     title(strcat('Frame n. ',num2str(i-1)));

% Panel first row
     subplot(4,4,1);
        title(strcat('Frame n. ',num2str(i-1)));
        imshow( cropped_image{1},[]);
        title(strcat('Phase Contrast'));
       
     subplot(4,4,2);
       imshow(imadjust(cropped_image{Parameters.fluoeval.fluorescence_channel }));
       % I tried to normalize the gfp,but there is not so much improvement.
       % Left as before.
%      gfp=(double(gfp-min_gfp))/double(max_gfp-min_gfp);
%      imshow(3*gfp);
       title(strcat('GFP'));
        
    subplot(4,4,3);
        imshow(imadjust(cropped_image{Parameters.fluoeval.nucleo_tag_channel}));
%       h2b=cropped_image{Parameters.fluoeval.nucleo_tag_channel};
%       h2b=(double(h2b-min_h2b))/double(max_h2b-min_h2b);
%       imshow(20*h2b);
        title(strcat('H2B'))
    
%     subplot(4,4,4);
% %       imshow(imadjust(cropped_image{4}))
%         imshow(imadjust(cropped_image{Parameters.sulpho}));
%         % sulpho=cropped_image{Parameters.sulpho};
%         % sulpho=(double(sulpho-min_sulpho))/double(max_sulpho-min_sulpho);
%         % imshow(4*sulpho);
% %       imread(fileNamesBCKGRND{4});
%         title(strcat('Sulforhodamine B'))
% 
% 
        %%
        % Panle second row
        
     subplot(4,4,5);
        imshow((Output.MASK{i}));
        title(strcat('Mask Frame'));
        
     subplot(4,4,6);
        daw1c=imread(fileNamesBCKGRND{1});
        imshow(imadjust(daw1c));
        title(strcat('DAW PH'));

     subplot(4,1,3);
     title(strcat('Frame n. ',num2str(i-1)));
%%    Choose the input plot
%     input=[ones(3*60,1)' zeros((n_timeframes-1-4)*60,1)']
      input=ones((Parameters.n_timeframes-1)*60,1)';
%     input=[ones(4*60,1)' zeros(24*60,1)' ones((Parameters.n_timeframes-1-4)*60,1)']
      plot(input,'color','[0.6350, 0.0780, 0.1840]','linewidth',1.5)
%%    Lines to define the experiment days 
%           xline(4*60,'color','r','linewidth',0.5,'LineStyle','--')
%           xline(4*60+24*60,'color','r','linewidth',0.5,'LineStyle','--')
%           xline(4*60+48*60,'color','r','linewidth',0.5,'LineStyle','--')
%           xline(4*60+24*3*60,'color','r','linewidth',0.5,'LineStyle','--')
            xlim([0 (n_timeframes-1)*60]);
            ylim([0 1])   
%           legend ('Input','Avarage fluo calibration phase ','Days off condition ''fontsize',6,'location','best')
            yticklabels({'-2i/Lif','+2i/Lif'})
            yticks([0 1])
            ylabel('Input','fontsize',10);
            title(strcat('Frame n. ',num2str(i-1)));  
            set(gca,'FontSize',10);
        
    subplot(4,1,4);
      cb=3; %value to moidify according to the frame of the calibration phase (usually 3hrs)
% normalizationp=(Output.Fluorescence_Channel_A_A(:,Parameters.fluoeval.fluorescence_channel)-min(Output.Fluorescence_Channel_A_A(:,Parameters.fluoeval.fluorescence_channel)))/(max(Output.Fluorescence_Channel_A_A(1:cb,Parameters.fluoeval.fluorescence_channel))- min(Output.Fluorescence_Channel_A_A(:,Parameters.fluoeval.fluorescence_channel)));

temp = Output.Fluorescence_Channel_A_A(:,Parameters.fluoeval.fluorescence_channel);
min_temp = min(temp);
max_temp = max(temp(1:cb));
normalizationp = (temp - min_temp) / (max_temp - min_temp);

      plot((0:i-1)*60,normalizationp(1:i),'linewidth',1.5,'color','[0.2 0.6 0]')
%     plot(Output.Nuclei_Normalized_Fluorescence_Channel_A_A(:,Parameters.fluoeval.fluorescence_channel))
%     average2ip=mean(normalizationp(1:cb,1));
%     average2i50p=average2ip*50/100;    
%     massimop=max(normalizationp(1:cb,1));
%     massimo50p=massimop*50/100;
      xlim([0 (Parameters.n_timeframes-1)*60]);
%     yline(average2ip,'color','[0.2 0.6 0]','linewidth',1,'LineStyle','--')
%     yline(average2i50p,'color','[0.2 0.6 0]','linewidth',0.5)
% %   yline(massimop,'color','[1 0.7 0.9]','linewidth',1)
% %   yline(massimo50p,'color','[1 0.7 0.9]','linewidth',0.5)   
      ylim([min(normalizationp) max(normalizationp)]);
%     legend ('GFP','<GFP> calibration phase','50% <GFP> calibration phase','fontsize',8,'location','best')
      legend ('GFP','fontsize',8,'location','best')
%     legend ('GFP/mCherry Normalized','Avarage fluo calibration phase ','50% meanfluo calibrationphase','Maximum calibration phase','50% maxfluo calibration phase','location','best','fontsize',12)
      ylabel('Fluo [a.u.]','fontsize',10);
      xlabel('Time [min]','fontsize',10);
      xlim([0 (Parameters.n_timeframes-1)*60]);
      set(gca,'FontSize',10)
      % hold off;

    % end
    pause(.2);
    frame = getframe(gcf);
    writeVideo(time_lapse,frame);
    pause(.1);
    writeVideo(time_lapse,getframe(fig1));

end
close(time_lapse); % close file for storing video data
end



