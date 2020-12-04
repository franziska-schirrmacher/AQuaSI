function E = calculateEnergyFunctionMC(f, g, grad, sigma, mu, lambda,...
    patchsize, quantile,guidance,weighting,dev,guided,Down,factor,weightColor, M)

% Calculate data fidelity term.
[~,~,numberFrames] = size(g);
[~,~,D] = size(guidance);
numberPixel = numel(f) / numberFrames;
r = zeros(numel(g),1);
% for d = 1:D
%     
%     switch weighting
%         case 'weighted'
%             %%Compute the weighted median
%             switch guided
%                 case 'intermediate'
%                     [LookUpRow, LookUpColumn] = main_aquasi_mex(image,image,patchsize,dev,quantile);
%                 case 'guidance'
%                     [LookUpRow, LookUpColumn] = main_aquasi_mex(image,guidance,patchsize,dev,quantile);
%             end
%         case 'unweighted'
%             %%Compute the median
%             [LookUpRow, LookUpColumn] = main_quasi_mex(image,guidance,patchsize,dev,quantile);
%     end
%     
%     M{d} = speye(numberPixel, numberPixel) - sparse(double(LookUpRow),double(LookUpColumn),1,numberPixel,numberPixel);
% end

E_tv = 0;
E_quasi = 0;

for i = 1:numberFrames
    start = (i-1)*numberPixel+1;
    ende = i*numberPixel;
    start2 = (i-1)*numberPixel/(factor*factor)+1;
    ende2 = i*numberPixel/(factor*factor);
    % Calculate residual error for Data Term.
    r(start2:ende2) = f(start:ende) - matrixToVector(g(:,:,i));
    
    % Calculate TV regularization term.
    for pos = 1:numel(grad)
        E_tv = E_tv + sum( abs(getGradient(f(start:ende),grad,pos)) );
    end
    % Calculate QuaSI regularization term.
    if D == numberFrames
        E_quasi = E_quasi + sum( abs(M{i} * f(start:ende)) );
    else
        E_quasi = E_quasi + sum( abs(M{1} * f(start:ende)) );
    end
end

% = repmat(Down*f(:),numberFrames,1) - g(:);
E_data = sum( huber(r, sigma) );

% Overall energy function.
E = E_data + mu*E_tv + lambda*E_quasi;
end

function h = huber(x, sigma)

h = 1/2 * x.^2;
h(abs(x) > sigma) = sigma * (abs( x(abs(x) > sigma) ) - 1/2*sigma);
end