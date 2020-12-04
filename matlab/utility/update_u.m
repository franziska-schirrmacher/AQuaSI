function [ u_vec ] = update_u( f_vec,M,bu_vec,lambda,alpha)

%Optimize u using soft thresholding (see equation (5) 
%please note, that M = Id - Q
u_vec = zeros(size(bu_vec));
if(alpha ~= 0)      
    %soft thresholding
    u_vec = shrink(( M * f_vec + bu_vec),lambda/alpha);    
end