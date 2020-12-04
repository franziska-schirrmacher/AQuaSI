% Model downsampling as matrix-vector product.
function D = getDownsampleMatrix(M, N, magFactor)

    numHRPixel = M*N;
    m = (M / magFactor);
    n = (N / magFactor);
    dindices = reshape(1:numHRPixel, M, N);
    dindices = dindices(1:magFactor:M, 1:magFactor:N);
    dindices = matrixToVector(dindices);

    D = sparse(1:(m*n), dindices', ones(size(dindices)), m*n, numHRPixel);