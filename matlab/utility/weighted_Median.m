function median = weighted_Median(patch,guide,center,dev)
%compute the weighted median for a given patch of the noisy image using
%the corresponding patch in the guidance image to compute the weights

intensityCenter = guide(center,center);
[R,C] = size(patch);
weights = zeros(R,C);



weightsTotal = 0;
for row = 1:R
    for col = 1:C        
        weights(row,col)  = exp(-(abs(double(guide(row,col)) ...
              - double(guide(center,center)))).^2/(2*(dev.^2)));
        weightsTotal = weightsTotal + weights(row,col);  
    end
end



% Reshape to vectors.
if ~isvector(weights)
    weights = weights(:);
end
if ~isvector(patch)
    patch = patch(:);
end

% Sort the data vector according to the given weight vector.
Image_and_weights_sorted = sortrows([patch weights], 1);
ImageSorted = Image_and_weights_sorted(:,1);
weightsSorted = Image_and_weights_sorted(:,2);

% Compute the cumulative sum of the weights.

% Select the weighted median according to the cumulative sum of the
% weights.
median = ImageSorted(1);
weights_curr  = 0;
for k = 1:length(weightsSorted)
    weights_curr = weights_curr  + weightsSorted(k);
    if  weights_curr >= 0.5*weightsTotal
        % We found the weighted median as the sum exceeds 0.5.
        median = ImageSorted(k);
        return;
    end
end
