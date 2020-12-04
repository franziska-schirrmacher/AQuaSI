close all;
clear all;
addpath(genpath('../algorithms/MCAQuaSI'));
addpath(genpath('../utility'));

%% specify the folders where data is stored and results are stored
dataDir = '../../data/bowls/';
resultDir = '../../results/bowls/';
if ~exist(resultDir, 'dir')
   mkdir(resultDir)
end

%% parameter for the experiment
maxOuterIter = 40;
maxInnerIter = 2;
maxCgIter = 20;


params.mu = 0.015;
params.lambda = 50;
params.alpha = 2000;
params.beta =0.3; 
params.maxOuterIter = 40;
params.maxInnerIter = 2;
params.maxCgIter = 50;
params.patchsize = 5;
params.dev = 0.1;


%% load the data
noisy =  im2double(imread(sprintf('%s/bowls1_noise.tiff',dataDir)));
NIR = im2double(imread(sprintf('%s/bowls1_ir.tiff',dataDir)));

rangeY = 1:size(NIR,1);
rangeX = 1:size(NIR,2);


noisy =  noisy(rangeY,rangeX,:);
NIR = NIR(rangeY,rangeX);

figure;imshow(noisy);title("input")
figure;imshow(NIR);title("NIR")


%% QuaSI
params.weighting = 'unweighted';

resultQuaSI = admmMC(noisy,params);

figure;imshow(resultQuaSI);title("QuaSI")

% Save results
save([resultDir sprintf('%s/QuaSI.mat',resultDir)],'resultQuaSI');
imwrite(resultQuaSI,[resultDir sprintf('%s/QuaSI.png',resultDir)]);

%% AQuaSI without guidance image using the input to compute Q
params.cases = 'guidance'; 
params.weighting = 'weighted';
params.guide = noisy; 

resultAQuaSIWOG_fixed = admmMC(noisy,params);    

figure;imshow(resultAQuaSIWOG_fixed);title("AQuaSI - Dynamic guidance, fixed f")

% Save results
save([resultDir sprintf('%s/AQuaSIWOG_fixed.mat',resultDir)],'resultAQuaSIWOG_fixed');
imwrite(resultAQuaSIWOG_fixed,[resultDir sprintf('%s/AQuaSIWOG_fixed.png',resultDir)]);


%% AQuaSI without guidance image using the intermediate results to compute Q
params.cases = 'intermediate'; 
params.weighting = 'weighted';

resultAQuaSIWOG = admmMC(noisy,params);    

figure;imshow(resultAQuaSIWOG);title("AQuaSI - Dynamic guidance, updated f")

% Save results
save([resultDir sprintf('%s/AQuaSIWOG.mat',resultDir)],'resultAQuaSIWOG');
imwrite(resultAQuaSIWOG,[resultDir sprintf('%s/AQuaSIWOG.png',resultDir)]);

%% AQuaSI multi-channel with guidance image available
params.guide = NIR;       
params.cases = 'guidance'; 
params.weighting = 'weighted';

resultAQuaSIWG = admmMC(noisy,params);     

figure;imshow(resultAQuaSIWG);title("AQuaSI - Static guidance")

% Save results
save([resultDir sprintf('%s/AQuaSIWG.mat',resultDir)],'resultAQuaSIWG');
imwrite(resultAQuaSIWG,[resultDir sprintf('%s/AQuaSIWG.png',resultDir)]);
    
    
    

