function [ result ] = cgCallMC( f_new_vec, wk, alpha, beta,...
    M, MT, L,Down,DownT,Dim,DGuide,weightColor,factor,params)
%use cg to solve the linear system of equations in equation (3.5) of the
%master's thesis


if any(strcmp('kernel',fieldnames(params)))
    ForwardFunc =  @(in_im) matrixToVector(imfilter(vectorToMatrix(in_im,params.R,params.C),rot90(params.kernel,2),'conv','same','replicate'));    
    BackwardFunc = @(in_im) matrixToVector(imfilter(vectorToMatrix(in_im,params.R,params.C),params.kernel,'conv','same','replicate'));    
else
    ForwardFunc = @(in_im) in_im;
    BackwardFunc = @(in_im) in_im;    
end

[totalPixel,~] = size(f_new_vec);
result = zeros(size(f_new_vec));
numPixel = totalPixel/ Dim;
for j = 1:Dim
    start = (j-1)*numPixel+1;
    ende = j*numPixel;
    start2 = (j-1)*numPixel/(factor*factor)+1;
    ende2 = j*numPixel/(factor*factor);
    
    if DGuide == Dim
        result(start:ende) = (2 *  BackwardFunc(DownT*(wk(start2:ende2) .*(Down*ForwardFunc(f_new_vec(start:ende)))))) ... % Data term  D^T*W^T*D.*f_new_vec
            + alpha*weightColor(j) *weightColor(j).* MT{j} * M{j} * f_new_vec(start:ende)...    % QuaSI prior
            + beta .* getLaplacian(f_new_vec(start:ende),L);
    else
        result(start:ende) = (2 * BackwardFunc(DownT*(wk(start2:ende2) .* (Down*ForwardFunc(f_new_vec(start:ende)))))) ... % Data term  D^T*W^T*D.*f_new_vec
            + alpha*weightColor(j) *weightColor(j).* MT{1} * M{1} * f_new_vec(start:ende)...    % QuaSI prior
            + beta .* getLaplacian(f_new_vec(start:ende),L);
    end
    
    % TV prior
end

end