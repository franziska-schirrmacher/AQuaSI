function [ result ] = cgCall( f_new_vec, wk, alpha, beta,...
    M, MT, L,D,DT,factor,params)
%use cg to solve the linear system of equations in equation (3.5) of the
%master's thesis


if any(strcmp('kernel',fieldnames(params)))
    ForwardFunc =  @(in_im) matrixToVector(imfilter(vectorToMatrix(in_im,params.R,params.C),rot90(params.kernel,2),'conv','same','replicate'));    
    BackwardFunc = @(in_im) matrixToVector(imfilter(vectorToMatrix(in_im,params.R,params.C),params.kernel,'conv','same','replicate'));    
else
    ForwardFunc = @(in_im) in_im;
    BackwardFunc = @(in_im) in_im;    
end


result = zeros(size(f_new_vec));
[numPixel,~] = size(f_new_vec);


start2 = 1;
ende2 = numPixel/(factor*factor);

result = (2 * BackwardFunc(DT*(wk(start2:ende2) .* (D*ForwardFunc(f_new_vec))))) ... % Data term  D^T*W^T*D.*f_new_vec
    + alpha .* MT * M * f_new_vec...    % QuaSI prior
    + beta .* getLaplacian(f_new_vec,L);                              % TV prior


end