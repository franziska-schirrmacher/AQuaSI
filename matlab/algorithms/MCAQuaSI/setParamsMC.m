function [params,f_new_vec,Down,DownT,E,grad,v_vec,bv_vec,u_vec,bu_vec,L,multiGuidance] = setParamsMC(g,params)

% regularization Parameter TV
if ~any(strcmp('mu',fieldnames(params)))
    params.mu = 0.05*3;
end
% regularization Parameter AQuaSI
if ~any(strcmp('lambda',fieldnames(params)))
    params.lambda = 13*3;
end
% Lagrangian multiplier AQuaSI
if ~any(strcmp('alpha',fieldnames(params)))
    if params.lambda == 0
        params.alpha = 0;
    else
        params.alpha = 1100*3;
    end
end
% Lagrangian multiplier TV
if ~any(strcmp('beta',fieldnames(params)))
    if params.mu == 0
         params.beta = 0;
    else
        params.beta = 7*3;
    end
end


% number of outer iteration
if ~any(strcmp('maxOuterIter',fieldnames(params)))
    params.maxOuterIter = 80;
end
% number of iterations with constant Q
if ~any(strcmp('maxInnerIter',fieldnames(params)))
    params.maxInnerIter = 1;
end

if ~any(strcmp('maxCgIter',fieldnames(params)))
    params.maxCgIter = 50;
end
% size of the kernel of the weighted median filter
if ~any(strcmp('patchsize',fieldnames(params)))
    params.patchsize = 5;
end
% standard deviation of the Gaussian weighting function for AQuaSI
if ~any(strcmp('dev',fieldnames(params)))
    params.dev = 0.05;
end

if ~any(strcmp('quantile',fieldnames(params)))
    params.quantile = 0.5;
end

% Decide whether to use AQuaSI (weighted) or QuaSI (unweighted) 
if ~any(strcmp('weighting',fieldnames(params)))
    params.weighting = 'weighted';
end

if ~(strcmp(params.weighting,'weighted') || strcmp(params.weighting,'unweighted'))
    error('weighting is either weighted (AQuaSI) or unweighted (QuaSI)');
end

% choose whether to use the guidance image or the intermediate results to
% compute the weights for the weighted quantile filter
if ~any(strcmp('cases',fieldnames(params)))
    params.cases = 'guidance';
end

if ~(strcmp(params.cases,'guidance') || strcmp(params.cases,'intermediate'))
    error('Case is either guidance or intermediate');
end

% super-resolution factor
if ~any(strcmp('factor',fieldnames(params)))
    params.factor = 1;
end

% if no guidance image is provided, the input is used as guidance image
if ~any(strcmp('guide',fieldnames(params)))    
    params.guide = imresize(g,params.factor,'bicubic');
end


% stop before max iterations is reached if there is no change between the
% previous and the current result
if ~any(strcmp('stopEarly',fieldnames(params)))
    params.stopEarly = true;
end
% choose your tolerance for early stopping
if ~any(strcmp('tol',fieldnames(params)))
    params.tol = 0.0001;
end



[Rn,Cn,Dim] = size(g);
[R,C,DGuid] = size(params.guide);

params.R = R;
params.C = C;

multiGuidance = false;
if DGuid == 3
    multiGuidance = true;
end


if (Dim ~= 3)
    error('Multi-channel approach -> input must be a RGB image!')
end


diff = C / Cn;
diff2 = R / Rn;
if(diff ~= params.factor || diff2 ~= params.factor)
    error('Downsampling factor must be an integer or the dimensions do not fit!')
end



% Initialization
if ~any(strcmp('initial',fieldnames(params)))    
    f_new_vec = matrixToVector(imresize(g,params.factor,'bicubic'));
else
    f_new_vec = matrixToVector(params.initial);
end

%get the downsampling matrix
if params.factor == 1
    Down = speye(R*C);
    DownT = Down;
else
    Down = getDownsampleMatrixExt(R,C,params.factor);
    DownT = (params.factor*params.factor)*Down';
end

%number of directions for the gradient
params.directions = 2;
% number of pixels per channel
numberPixel = R*C;

%define the gradient directions and compute its transposes
gradX = getGradientMatrix(numberPixel,R);
gradY = getGradientMatrix(numberPixel,1);
grad = {gradX, gradY};

clearvars gradX;
clearvars gradY;
E = [];

%iniitalize the auxiliary variable and the bregman variable related to the
%total variation term
v_vec = cell(1,params.directions);
bv_vec = cell(1,params.directions);

for pos = 1:params.directions
    v_vec{pos} = matrixToVector(zeros(R,C*Dim));
    bv_vec{pos} = matrixToVector(zeros(R,C*Dim));
end

%iniitalize the auxiliary variable and the bregman variable related to the
%quasi prior
u_vec = zeros(R*C*Dim,1);
bu_vec = zeros(R*C*Dim,1);


%compute Laplacian
[Rgrad, Cgrad] = size(grad{1});
L = sparse(Rgrad,Cgrad);
for pos  = 1:length(grad)
    L = L + grad{pos}'*grad{pos};
end




end

