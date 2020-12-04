close all;
clear all;
addpath(genpath('../algorithms/MCAQuaSI'));
addpath(genpath('../utility'));

%% specify the folders where data is stored and results are stored
dataDir = '../../data/ARRI/';
resultDir = '../../results/ARRI/';

if ~exist(resultDir, 'dir')
   mkdir(resultDir)
end

dataset = dir(dataDir);

%% parameter for the experiment
maxOuterIter = 40;
maxInnerIter = 2;
maxCgIter = 20;


frames = 1;
params.mu = 0.05*3;
params.lambda = 13*3;
params.alpha = 1100*3;
params.beta =7*3;
params.updateMultiplier = false;   
params.maxOuterIter = 40;
params.maxInnerIter = 2;
params.maxCgIter = 50;



for index = 1:numel(dataset)
    
    if strcmp(dataset(index).name, '.') || strcmp(dataset(index).name, '..')
        continue
    end
    %% load the data from the data folder
    RGB_name = dir([dataDir, '/', dataset(index).name, '/RGB','*']);
    RGB = im2double(imread([dataDir '/' dataset(index).name '/' RGB_name.name]));
    NIR_name = dir([dataDir, '/', dataset(index).name, '/NIR','*']);
    NIR = im2double(imread([dataDir '/' dataset(index).name '/' NIR_name.name]));
    rangeY = 1:size(NIR,1);
    rangeX = 1:size(NIR,2);
    name = RGB_name.name;    

    noisy =  RGB(rangeY,rangeX,:);    
    noisy = imnoise(noisy,'speckle',0.04);
    RGB =  RGB(rangeY,rangeX,:);
    NIR = NIR(rangeY,rangeX);
    
    %% print the name of the current file
    parseName = strsplit(RGB_name.name,'.tif');
    parseName = parseName{1};
    parseName = strsplit(parseName,'RGB_');
    parseName = parseName{2};
    parseName = strsplit(parseName,'.');
    parseName = parseName{1};
    disp(parseName)



    %% AQuaSI without guidance
    params.cases = 'intermediate'; 
    params.patchsize = 5;
    params.dev = 0.02;
    resultAQuaSIWOG = admmMC(noisy,params);     
        
    % Save results
    save([resultDir sprintf('%s/%s_MCWOG.mat',resultDir,parseName)],'resultAQuaSIWOG');
    imwrite(resultAQuaSIWOG,[resultDir sprintf('%s/%s_MCWOG.png',resultDir,parseName)]);

    %% AQuaSI multi-channel with guidance
    params.guide = NIR;       
    params.cases = 'guidance'; 
    params.patchsize = 9;
    params.dev = 0.006;

    resultAQuaSIWG = admmMC(noisy,params);     
        
    % Save results
    save([resultDir sprintf('%s/%s_MCWG.mat',resultDir,parseName)],'resultAQuaSIWG');
    imwrite(resultAQuaSIWG,[resultDir sprintf('%s/%s_MCWG.png',resultDir,parseName)]);
    
    
    

end
