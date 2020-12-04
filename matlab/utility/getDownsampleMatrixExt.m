function D = getDownsampleMatrixExt(M, N, magFactor)

	avg = 1 / (magFactor*magFactor);

	numHRPixel = M*N;
    m = (M / magFactor);
    n = (N / magFactor);
    dindices = reshape(1:numHRPixel, M, N);
    dindices = dindices(1:magFactor:M, 1:magFactor:N);
    dindices = matrixToVector(dindices);
	
    row = [];
    
    for j = 1:n
        rowTmp = [];
        for i = 1:m
            rowTmp = [rowTmp, i+((j-1)*m)];            
        end
        row = [row, repmat(rowTmp,1,magFactor)];
    end
    
   row = repmat(row,1, magFactor);
%     row = [];
%     for i = 1:m*n
%          row = [row, i];
%     end
%     row = repmat(row,1,magFactor*size(row,1));
%     row = repmat(row,1,magFactor);
    
	column = 1:magFactor:numHRPixel;
	
    for i = 2:magFactor	
		column = [column,i:magFactor:(numHRPixel)];
    end
   
    D = sparse(row, column, ones(1,(m*n)*(magFactor*magFactor)).*avg, m*n, numHRPixel);
end