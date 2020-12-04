function [f_new_vec, sigma] = update_f( f_old_vec,g,alpha,beta,v_vec,bv_vec,bu_vec,u_vec,...
    M,MT,grad,L,tol,maxCgIter,Down,DownT,factor,params)

% Assemble right hand side of the equation system.
% Part for TV denoising and the QuaSI prior

if any(strcmp('kernel',fieldnames(params)))
    ForwardFunc =  @(in_im) matrixToVector(imfilter(vectorToMatrix(in_im,params.R,params.C),rot90(params.kernel,2),'conv','same','replicate'));    
    BackwardFunc = @(in_im) matrixToVector(imfilter(vectorToMatrix(in_im,params.R,params.C),params.kernel,'conv','same','replicate'));    
else
    ForwardFunc = @(in_im) in_im;
    BackwardFunc = @(in_im)in_im;    
end


[~,directions] = size(v_vec);
[R,C] = size(v_vec{1});
TVpart =zeros(R*C,1);


% TV part
if beta > 0
    for pos = 1: directions
        TVpart = TVpart...
            + getGradient(bv_vec{pos} - v_vec{pos}...
            ,grad,pos, 'backward');
    end
end
% Determine residual error.
r = vectorToMatrix(Down*ForwardFunc(f_old_vec), size(g,1), size(g,2)) - g;
if (all(matrixToVector(r) == 0))
    % No adaptive noise estimation.
    sigma = Inf;
    w = ones(size(r));
else
    % Adpative noise estimation.
    res = r;
    sigma = 1.4826 * mad(res(:), 1);
    
    if(sigma == 0)
        w = ones(size(r));
    else
        res(res == 0) = (1.345*sigma);
        r = res;
        w = (1.345*sigma) ./ abs(r);
    end
    
    % Compute confidence weights according to the Huber loss.
    % Set threshold to 95-percent efficiency
    
end
w(abs(r) <= 1.345*sigma) = 1;

b = (2*BackwardFunc(DownT*matrixToVector(w .* g))) - beta *TVpart ...
    - alpha * MT * (bu_vec - u_vec);


% Assemble matrix (implemented via imfilter)
[Af] = @(f_new_vec)cgCall( f_new_vec, matrixToVector(w)...
    , alpha, beta, M, MT, L,Down,DownT,factor,params);



% Solve A * f_new = b
[f_new_vec,~,~,~,~] = cgs(Af, b, ...
    tol,maxCgIter, [], [],...
    f_old_vec);


% Truncate intensity values to [0, 1].
f_new_vec(f_new_vec < 0) = 0;
f_new_vec(f_new_vec > 1) = 1;

end