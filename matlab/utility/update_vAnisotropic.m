function [v ] = update_vAnisotropic(  gradient, bv,beta,mu  )
% Update the auxiliary variable v according to equation (6)
v = zeros(size(gradient));

if(beta ~= 0)    
   v = shrink((gradient+bv),mu/(beta));
end

end

