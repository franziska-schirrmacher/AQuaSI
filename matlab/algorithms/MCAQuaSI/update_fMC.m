function [f_new_vec, sigma] = update_fMC( f_old_vec,g,alpha,beta,v_vec,bv_vec,bu_vec,u_vec,...
    M,MT,grad,L,tol,maxCgIter,Down,DownT,Dim,DGuide,weightColor,factor,params)%,gradient,mu,lambda,patchsize,quantile,goldstandard)

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
b = zeros(size(bu_vec));
numPixel = R*C/Dim ;




for j = 1:Dim
    start = (j-1)*numPixel+1;
    ende = j*numPixel;
    % TV part
    if beta > 0
        for pos = 1: directions
            TVpart(start:ende) = TVpart(start:ende)...
                + getGradient(bv_vec{pos}(start:ende) - v_vec{pos}(start:ende)...
                ,grad,pos, 'backward');
        end
    end
    % Determine residual error.
    r(:,:,j) = vectorToMatrix(Down*ForwardFunc(f_old_vec(start:ende)), size(g,1), size(g,2)) - g(:,:,j);
    if (all(matrixToVector(r(:,:,j)) == 0))
        % No adaptive noise estimation.
        sigma = Inf;
        w(:,:,j) = ones(size(r(:,:,j)));
    else
        % Adpative noise estimation.
        res = r(:,:,j);
        sigma = 1.4826 * mad(res(:), 1);
        
        if(sigma == 0)
            w(:,:,j) = ones(size(r(:,:,j)));
        else
            res(res == 0) = (1.345*sigma);
            r(:,:,j) = res;
            w(:,:,j) = (1.345*sigma) ./ abs(r(:,:,j));
        end
        
        % Compute confidence weights according to the Huber loss.
        % Set threshold to 95-percent efficiency
        
    end
    w(abs(r(:,:,j)) <= 1.345*sigma) = 1;
    if DGuide  == Dim
        b(start:ende) = (2*BackwardFunc(DownT*(matrixToVector(w(:,:,j) .* g(:,:,j))))) - beta *TVpart(start:ende) ...
            - alpha*weightColor(j) * MT{j} * (bu_vec(start:ende) - u_vec(start:ende));
    else
        b(start:ende) = (2*BackwardFunc(DownT*(matrixToVector(w(:,:,j) .* g(:,:,j))))) - beta *TVpart(start:ende) ...
            - alpha*weightColor(j) * MT{1} * (bu_vec(start:ende) - u_vec(start:ende));
    end
    
    
    
end


% Assemble matrix (implemented via imfilter)
[Af] = @(f_new_vec)cgCallMC( f_new_vec, matrixToVector(w)...
    , alpha, beta, M, MT, L,Down,DownT,Dim,DGuide,weightColor,factor,params);



% Solve A * f_new = b
[f_new_vec,~,~,~,~] = cgs(Af, b, ...
    tol,maxCgIter, [], [],...
    f_old_vec);


% Truncate intensity values to [0, 1].
f_new_vec(f_new_vec < 0) = 0;
f_new_vec(f_new_vec > 1) = 1;

end