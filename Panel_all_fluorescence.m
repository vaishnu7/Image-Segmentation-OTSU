%% Panel figures showing fluorescence   

%% General info on the variables name (for more info see .mat fluorescence_evaluation)

%% Output.Absolute_fluorescence_AVG (absolute fluo as avg of pixels)
% AVG_measure_from_MaskedFrame_Channel=mean(Masked_Frame_from_Channel{i}(Masked_Frame_from_Channel{i}~=0))  

%% Output.Fluorescence_Channel_A_A (avg of masked fluo frame - avg background frame) 
%   AVG_measure_from_MaskedFrame_Channel(i) - AVG_measure_from_BCKGRND_Frame_Channel(i)

% Or (in case you make computation with the antimask option)

%   AVG_measure_from_MaskedFrame_Channel(i) -  AVG_measure_from_antimask_BCKGRND_Channel(i)

%% Output.Average_Background_fluorescence
%   AVG_measure_from_BCKGRND_Frame_Channel(i) = mean(mean(BCKGRND_Frame_from_Channel{i}));

% or (in case you make computation with the antimask option)

%   AVG_measure_from_antimask_BCKGRND_Channel(i) = mean(BCKGRND_from_antimask{i}(BCKGRND_from_antimask{i}~=0))

%% Output.Nuclei_Normalized_Fluorescence_Channel_A_A
% Fluorescence_Channel_A_A(i)/Fluorescence_Channel_A_A(Input_data.nucleo_tag_channel)


%% General variables allocation
  cf=2; %GFP (CHANNEL GFP FLUO)
  cb=3; %CALIBRATION PHASE (NUMBER OF FRAME)
  cn=4; %H2B (CHANNEL NUCLEO TAG)
  input=ones((Parameters.n_timeframes)*60,1);   
%%
figure

%% input
subplot(2,2,1)
plot(input,'color','[0.6350, 0.0780, 0.1840]','linewidth',1)
xlim([0 (Parameters.n_timeframes-1)*60]);
yticks([0 1])
% yticklabels({'Ndiff','2i'})
yticklabels({'-2i','2i'})
% yticklabels({'-(2i+LIF)','(2i+LIF)'})
ylabel('Input','fontsize',12);
% title(' 2i 100920  pos 15')
ylim([-0.1 1.1])
set(gca,'FontSize',12)

%% fluorescence background
subplot(2,2,2)
average2i=mean(outputArg1.Average_Background_fluorescence (1:cb,cf));
average2i50=average2i*50/100;
plot(0:60:(Parameters.n_timeframes-1)*60,outputArg1.Average_Background_fluorescence(:,cf),'color','g','linewidth',1);
hold on
% plot(0:60:(Parameters.n_timeframes-1)*60,outputArg1.Average_Background_fluorescence(:,cn),'color','k','linewidth',1);
% legend ('Average Background fluorescence GFP','Average Background fluorescence IRP','location','best','fontsize',6)
legend ('Average Background fluorescence GFP','location','best','fontsize',6);
ylabel('Fluo [A.u.]','fontsize',12);
% xlabel('Time [min]','fontsize',12);
xlim([0 (Parameters.n_timeframes-1)*60]);
set(gca,'FontSize',12)


%% fluorescence GFP and IRFP
subplot(2,2,3);
% average2i=mean(outputArg1.Absolute_Fluorescence_AVG (1:cb,cf));
% average2i50=average2i*50/100;
plot(0:60:(Parameters.n_timeframes-1)*60,outputArg1.Absolute_Fluorescence_AVG (:,cf),'color','g','linewidth',1);
hold on;
plot(0:60:(Parameters.n_timeframes-1)*60,outputArg1.Absolute_Fluorescence_AVG (:,cn),'color','r','linewidth',1);
xlim([0 (Parameters.n_timeframes-1)*60]);
legend ('Absolute Fluorescence AVG  GFP','Absolute Fluorescence AVG IRFP','location','best','fontsize',6);
ylabel('Fluo [A.u.]');
xlabel('Time [min]','fontsize',12);
set(gca,'FontSize',12)

%% fluorescenza della gfp normalizzata per l'h2b
% subplot(3,2,4);
% average2i=mean(outputArg1.Nuclei_Normalized_Fluorescence_Channel_A_A (1:cb,cf));
% average2i50=average2i*50/100;
% plot(0:60:(Parameters.n_timeframes-1)*60,outputArg1.Nuclei_Normalized_Fluorescence_Channel_A_A(:,cf),'color','[0.9290, 0.6940, 0.1250]','linewidth',1);
% hold on
% xlim([0 (Parameters.n_timeframes-1)*60]);
% ylim([0 max(outputArg1.Nuclei_Normalized_Fluorescence_Channel_A_A (:,cf))+0.5]);
% yline(average2i,'color','[0.9290, 0.6940, 0.1250]','linewidth',1,'LineStyle','--')
% yline(average2i50,'color','[0.9290, 0.6940, 0.1250]','linewidth',0.5)
% legend ('Normalized channels','fontsize',6)
% ylabel('Fluo [A.u.]','fontsize',12);
% set(gca,'FontSize',12)

%% fluorescence channel A_A - GFP
subplot(2,2,4);
average2i = mean(outputArg1.Fluorescence_Channel_A_A  (1:cb,cf));
% average2i50=average2i*50/100;
% average2ir=mean(outputArg1.Fluorescence_Channel_A_A  (1:cb,cn));
plot(0:60:(Parameters.n_timeframes-1)*60,outputArg1.Fluorescence_Channel_A_A (:,cf),'color','[0.2 0.4 0]','linewidth',1);
%hold on
% plot(0:60:(Parameters.n_timeframes-1)*60,outputArg1.Fluorescence_Channel_A_A (:,cn),'color','r','linewidth',1);
yline(average2i,'color','[0.2 0.4 0]','linewidth',1,'LineStyle','--')
% yline(average2ir,'color','r','linewidth',1,'LineStyle','--')
% yline(average2i50,'color','k','linewidth',1,'LineStyle','--')
% legend ('Fluorescence Channel AA GFP','Fluorescence Channel AA IRP','Mean GFP','Mean IRFP','50%','fontsize',4,'location','best')
legend ('Fluorescence Channel AA GFP','Mean GFP','fontsize',4,'location','best')
ylabel('Fluo [A.u.]','fontsize',12);
xlabel('Time [min]','fontsize',12);
xlim([0 (Parameters.n_timeframes-1)*60]);
set(gca,'FontSize',12)
