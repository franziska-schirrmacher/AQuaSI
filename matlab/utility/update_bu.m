function [ bu_vec] = update_bu( bu_old_vec,f_vec,u_vec,M)
% update the bregman variable bu according to equation (7) 

bu_vec = bu_old_vec + M*f_vec - u_vec;

end

