function g_vec = getLaplacian(f_vec,L)
%compute the laplacian of the denoised image f (in vector form)

g_vec = L*f_vec;

