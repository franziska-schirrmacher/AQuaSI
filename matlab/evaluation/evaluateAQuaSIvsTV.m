close all;
clear all;

% path to the ARRI Dataset
dataDir = '../../data/ARRI/';
% path to the folder where the results are stored
resultDir = '../../results/AQuaSIvsTV';

if ~exist(resultDir, 'dir')
   mkdir(resultDir)
end

% get the dataset
dataset = dir(dataDir);

addpath(genpath('../algorithms/SCAQuaSI'));
addpath(genpath('../utility'));


% choose whether AQuaSI and/or TV are evaluated
TV = 1;
AQuaSI = 1;

% use AQuaSI without the a guidance image
cases = 'intermediate';
mu = 0.4;
lambda = 23;
alpha = 2000;
beta = 30;

iter = 1;
num = 2;
multiA = logspace(-2,3,num);
multiT = logspace(-3,1,num);
pA = zeros(num,1,4);
l = zeros(num,1);
pT = zeros(num,1);
m = zeros(num,1);

for index = 3%1:numel(dataset)
    
    if strcmp(dataset(index).name, '.') || strcmp(dataset(index).name, '..')
        continue
    end
    
    RGB_name = dir([dataDir, '/', dataset(index).name, '/RGB','*']);
    RGB = im2double(imread([dataDir '/' dataset(index).name '/' RGB_name.name]));
    NIR_name = dir([dataDir, '/', dataset(index).name, '/NIR','*']);
    NIR = im2double(imread([dataDir '/' dataset(index).name '/' NIR_name.name]));
  
    rangeY = 400:800;
    rangeX =  400:800;
    name = RGB_name.name;   

    
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
    
    initial = noisy;
    
    [R,C] = size(noisy);
    
    resultTV = zeros(R,C,num);
    muA = [0,0.0001,0.001,0.01, 0.1];
    betaA = [0,0.0075,0.075,0.75, 7.5];
    for j = 1:5
        resultAQuaSI = zeros(R,C,num);
        for i = 1:num
            params.maxOuterIter = 40;
            params.maxInnterIter = 2;
            params.maxCgIter = 10;
            % AQuaSI
            if AQuaSI
                
                params.lambda = lambda*multiA(i);
                params.alpha = alpha*multiA(i);
                params.mu = muA(j);
                params.beta = betaA(j);
                params.patchsize = 3;
                params.dev = 0.1;
                params.cases = cases;
                params.updateMultiplier = false;
                resultAQuaSI(:,:,i) = admmSC(noisy,params);
                
                
                pA(i,iter,j) = psnr(resultAQuaSI(:,:,i),RGB);
                l(i,iter) = lambda*multiA(i);                
            end            
            
            if TV && j == 1
                params.lambda = 0;
                params.beta = beta*multiT(i);
                params.mu = mu*multiT(i);
                params.alpha = 0;
                params.updateMultiplier = false;
                resultTV(:,:,i) = admmSC(noisy,params);
                
                pT(i,iter) = psnr(resultTV(:,:,i),RGB);
                m(i,iter) = mu*multiT(i);
            end
            
            
        end
        if AQuaSI
            save(sprintf('%s/resultAQuaSImu%d.mat',resultDir,j), 'resultAQuaSI');
        end
    end
    if TV
        save(sprintf('%s/resultTV.mat',resultDir), 'resultTV');
    end
    iter = iter + 1;
end




