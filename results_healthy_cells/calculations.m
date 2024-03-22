close all;
clear;
clc;

%% parameters to define

cf = 2; % GFP channel
cb = 3; % calibration phase
n_timeframes = 43 ; % time frame size
% average = zeros(n_timeframes,8); % to create an array for storing the average values
store = zeros(n_timeframes,8); % to create an array for storing the extracted Fluorescence_Channel_A_A

%% for position 3
n=1;
data = load("Results---3.mat");
store(:,n) = data.outputArg1.Fluorescence_Channel_A_A(:,2);
% average(n) = mean(data.outputArg1.Fluorescence_Channel_A_A(1:cb,cf));

%% for position 6
n=2;
data = load("Results---6.mat");
store(:,n) = data.outputArg1.Fluorescence_Channel_A_A(:,2);
% average(n) = mean(data.outputArg1.Fluorescence_Channel_A_A(1:cb,cf));

%% for position 7
n=3;
data = load("Results---7.mat");
store(:,n) = data.outputArg1.Fluorescence_Channel_A_A(:,2);
% average(n) = mean(data.outputArg1.Fluorescence_Channel_A_A(1:cb,cf));

n=n+1;

%% for position 10 to 14
for t = 10:1:14
    % name_position = strcat('Folder ', pad(num2str(t)), '/Results.mat');
    % if you want to redirect to the folder and not just the .mat file

    name_position = strcat('Results','---', pad(num2str(t)), '.mat'); % not convenient for large set of results
    data = load(name_position);
    store(:,n) = data.outputArg1.Fluorescence_Channel_A_A(:,2);
    % average(n) = mean(data.outputArg1.Fluorescence_Channel_A_A(1:cb,cf));
    n=n+1;
end
% calculate mean for all the rows (time) by all the healthy cells
% (number of healthy chambers)
n = 1;
for i=1:n_timeframes
average(i) = mean(store(i,8));
end
% standard deviation for each row
S = std (store,0,2);
SE = S/sqrt(length(store)); % standard error of the mean

%%
% figure
% 
% 
plot(average,'Linewidth', 1.5);
hold on
plot(S, 'b--o', 'Linewidth', 1.5);
hold on
errorbar(S,SE,'r--*', 'Linewidth', 1.2)
ylabel('Avg Fluo [A.u.]','fontsize',10);
% xticklabels({'n = 3','n = 6','n = 7','n = 10','n = 11','n = 12','n = 13','n = 14'})
% xlabel('Healthy Chambers','fontsize',10);
xlabel('Time [h]','FontSize',10)
legend ('mean','standard deviation with error bar');
grid on