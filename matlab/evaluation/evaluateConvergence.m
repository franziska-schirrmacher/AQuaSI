% IMPORTANT
% Computing the energy in every iteration corrupts the timing computation.
% The timings from the paper are provided

close all;
clear all;
addpath(genpath('../algorithms/SCAQuaSI'));
addpath(genpath('../utility'));

%% specify the folders where data is stored and results are stored
dataDir = '../../data/ARRI/';
resultDir = '../../results/convergence';

if ~exist(resultDir, 'dir')
   mkdir(resultDir)
end

dataset = dir(dataDir);


% Variables to store the value of the energy function
E1_all = [];
E2_all = [];

for index = 1:numel(dataset)
    
    if strcmp(dataset(index).name, '.') || strcmp(dataset(index).name, '..')
        continue
    end
    %% load the data from the data folder
    RGB_name = dir([dataDir, '/', dataset(index).name, '/RGB','*']);
    RGB = im2double(imread([dataDir '/' dataset(index).name '/' RGB_name.name]));
    NIR_name = dir([dataDir, '/', dataset(index).name, '/NIR','*']);
    NIR = im2double(imread([dataDir '/' dataset(index).name '/' NIR_name.name]));
    %% if this takes too long, you might want to work on smaller images patches
    rangeY = 1:size(NIR,1);
    rangeX = 1:size(NIR,2);
    name = RGB_name.name;
    
    % For this experiment, we work on grayscale data
    noisy =  RGB(rangeY,rangeX,:);    
    noisy = imnoise(noisy,'speckle',0.04);
    noisy = rgb2gray(noisy);
    RGB =  RGB(rangeY,rangeX,:);
    NIR = NIR(rangeY,rangeX);    
    RGB = rgb2gray(RGB);
   
    parseName = strsplit(RGB_name.name,'.tif');
    parseName = parseName{1};
    parseName = strsplit(parseName,'RGB_');
    parseName = parseName{2};
    parseName = strsplit(parseName,'.');
    parseName = parseName{1};
    disp(parseName)
    
    params.updateMultiplier = false;   
    params.maxCgIter = 50;
    params.patchsize = 9;
    params.dev = 0.006;
    params.guide = NIR;
    
    %% Compute Q in every iteration
    params.maxOuterIter = 80;
    params.maxInnerIter = 1;
    [resultAQuaSI,E1] = admmSC(noisy,params);
    E1_all = [E1_all; E1];
    
   
    %% Compute Q in every second iteration    
    params.maxOuterIter = 40;
    params.maxInnerIter = 2;
    
    [resultAQuaSI2,E2] = admmSC(noisy,params);
    E2_all = [E2_all; E2];
    
end



save(sprintf('%s/E1_all.mat',resultDir),'E1_all');
save(sprintf('%s/E2_all.mat',resultDir),'E2_all');

load(sprintf('%s/T1.mat',resultDir));
load(sprintf('%s/T2.mat',resultDir));

E1 = mean(E1_all,1);
E2 = mean(E2_all,1);

figure;
semilogy(T1,E1)
hold on
semilogy(T2,E2)
grid on;
legend('Update every iterations', 'Update every second iteration')




