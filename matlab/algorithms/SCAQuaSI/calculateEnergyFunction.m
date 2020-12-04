function E = calculateEnergyFunction(f, g, grad, sigma, mu, lambda,...
    patchsize, quantile,guidance,weighting,dev,guided,Down,factor)

% Calculate data fidelity term.
[R,C,Dim] = size(g);
numberPixel = numel(f) / Dim;
r = zeros(numel(g),1);
switch weighting
    case 'weighted'
        %%Compute the weighted median of the neighbourhoof
        switch guided
            
            case 'intermediate'
                [LookUpRow, LookUpColumn] = main_aquasi_mex(vectorToMatrix(f,R,C,Dim),vectorToMatrix(f,R,C,Dim),patchsize,dev,quantile);
            case 'guidance'
                [LookUpRow, LookUpColumn] = main_aquasi_mex(vectorToMatrix(f,R,C,Dim),guidance,patchsize,dev,quantile);
        end
        
    case 'unweighted'
        %%Compute the median of the neighbourhood
        [LookUpRow, LookUpColumn] = main_quasi_mex(vectorToMatrix(f,R,C,Dim),guidance,patchsize,dev,quantile);
end

M = speye(numberPixel, numberPixel) - sparse(double(LookUpRow),double(LookUpColumn),1,numberPixel,numberPixel);
E_tv = 0;
E_quasi = 0;

for i = 1:Dim
    start = (i-1)*numberPixel+1;
    ende = i*numberPixel;
    start2 = (i-1)*numberPixel/(factor*factor)+1;
    ende2 = i*numberPixel/(factor*factor);
    % Calculate residual error for Data Term.
    r(start2:ende2) = Down*f(start:ende) - matrixToVector(g(:,:,i));
    
    % Calculate TV regularization term.
    for pos = 1:numel(grad)
        E_tv = E_tv + sum( abs(getGradient(f(start:ende),grad,pos)) );
    end
    % Calculate QuaSI regularization term.
    E_quasi = E_quasi + sum( abs(M * f(start:ende)) );
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