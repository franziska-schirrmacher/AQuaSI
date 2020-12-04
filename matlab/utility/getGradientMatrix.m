function Q = getGradientMatrix(numPixel, delta)
% creates a matrix which yields the gradient of the image (image must be in
% vector form). The variable delta determines whether the gradient is
% computed in x-direction (delta = 1) or in y-direction (delta = #columns of the image)

    % Assemble pixel positions in linearized order.
    pixelPos = 1:numPixel;
    
    % Get positions (circular) shifted by a certain number of pixels. 
    pixelPos_delta = mod(pixelPos + delta, numPixel);
    pixelPos_delta(pixelPos_delta == 0) = numPixel;
    
    % Sparse matrix with -1 on a certain diagonal.
    V = sparse(pixelPos', pixelPos_delta, -ones(1,numPixel), numPixel, numPixel);
    % Sparse identity matrix.
    I = speye(numPixel, numPixel);
    % Get overall gradient matrix.
    Q = I + V;