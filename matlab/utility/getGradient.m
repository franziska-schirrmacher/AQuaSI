function [g] = getGradient(f_vec,grad,pos, dir)
%functions computes the gradient of the vectorized image f 
%the variable grad stores all the gradient matrices computed in admm.m
% the variable gradT stores the transposed matrices
% pos is the current gradient direction 
% dir determines whether grad or gradT is used

if nargin < 5
    dir = 'forward';
end

if strcmp(dir, 'forward')  
        
   g = grad{pos}*f_vec;   
    
elseif strcmp(dir, 'backward')      
        
    g = grad{pos}'*f_vec;  
    
end
