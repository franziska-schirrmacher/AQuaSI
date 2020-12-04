% Model blurring related to camera PSF as matrix-vector product.
function K = getBlurMatrix(M, N, blur_size, psfWidth)

    numHRPixel = M*N;
    sz_blur = blur_size * blur_size;
    ref_ind = floor(sz_blur / 2) + 1;

    % Assemble Gaussian blur kernel.
    H = fspecial('gaussian', [blur_size, blur_size], psfWidth);
    h = H(:);
    h = h';

    % Indicies
    rep = (1:blur_size)';
    x = rep(:, ones(blur_size,1)); 
    x = x(:)';
    rep = rep';
    y = rep(ones(blur_size,1),:); 
    y = y(:)';
    ind_h = sub2ind([M, N], x, y);
    ind_h_all = ind_h(ones(M*N,1), :)'; 
    ind_h_all = ind_h_all(:);

    % All pixels
    all = 1:M*N;
    all_ind = all(ones(sz_blur,1),:); 
    all_ind = all_ind(:);

    % Offset
    offset = all_ind - ind_h(ref_ind);
    new_ind = ind_h_all + offset;
    h_mat = mod(find(ind_h_all),sz_blur);
    h_mat(h_mat == 0) = sz_blur;
    h_mat = h(h_mat); 
    h_mat = h_mat(:);

    % Elimination of false indices
    false_ind = find(new_ind < 1);
    false_ind = cat(1, false_ind, find(new_ind > numHRPixel));
    all_ind(false_ind) = [];
    new_ind(false_ind) = [];
    h_mat(false_ind) = [];

    % Assemble sparse matrix
    K = sparse(all_ind, new_ind, h_mat, numHRPixel, numHRPixel);