function [ f,E] = admmMC(g,params)
%store the size of the image

[params,f_new_vec,Down,DownT,E,grad,...
    v_vec,bv_vec,u_vec,bu_vec,L,multiGuidance] = setParamsMC(g,params);

[~,~,Dim] = size(g);
[R,C,DGuid] = size(params.guide);
numberPixel = R*C;

if multiGuidance
    weightColor = [1,1,1];
else
    weightColor = [0.299 , 0.587 , 0.114 ];
end

for outerIter = 1:params.maxOuterIter
    %figure;imshow(vectorToMatrix(f_new_vec,R,C,Dim))
    if params.lambda > 0
        % Compute look up table for quantile regularization (M = Id - Q)
        
        image = vectorToMatrix(f_new_vec,R,C,Dim);
        if ~multiGuidance
            for i = 1:Dim
                image(:,:,i) =  weightColor(i)*image(:,:,i);
            end
            image = double(sum(image,3));
        end
        
        for d = 1:DGuid
            
            switch params.weighting
                case 'weighted'
                    %%Compute the weighted median
                    switch params.cases
                        case 'intermediate'
                            [LookUpRow, LookUpColumn] = main_aquasi_mex(...
                                image(:,:,d),image(:,:,d),...
                                params.patchsize,params.dev,params.quantile);
%                             [LookUpRow, LookUpColumn] = compute_Q(...
%                                 image(:,:,d) ,params.patchsize,params.quantile,...
%                                 image(:,:,d),params.weighting,params.dev, ...
%                                 params.cases,weightColor);
                        case 'guidance'
                            [LookUpRow, LookUpColumn] = main_aquasi_mex(...
                                image(:,:,d),params.guide(:,:,d),...
                                params.patchsize,params.dev,params.quantile);                            
%                             [LookUpRow, LookUpColumn] = compute_Q(...
%                                 image(:,:,d) ,params.patchsize,params.quantile,...
%                                 params.guide(:,:,d),params.weighting,params.dev, ...
%                                 params.cases,weightColor);                    
                    end
                case 'unweighted'
                    %%Compute the median
                    [LookUpRow, LookUpColumn] = main_quasi_mex(...
                        image(:,:,d),params.guide(:,:,d),...
                        params.patchsize,params.dev,params.quantile);
%                     [LookUpRow, LookUpColumn] = compute_Q(...
%                         image(:,:,d) ,params.patchsize,params.quantile,...
%                         image(:,:,d),params.weighting,params.dev, ...
%                         params.cases,weightColor);
            end
            
            M{d} = speye(numberPixel, numberPixel) - ...
                sparse(double(LookUpRow),double(LookUpColumn),1,numberPixel,numberPixel);
            MT{d} = speye(numberPixel, numberPixel) - ...
                sparse(double(LookUpColumn),double(LookUpRow),1,numberPixel,numberPixel);
        end
    else
        if ~multiGuidance
            M{1} = speye(numberPixel, numberPixel);
            MT{1} = speye(numberPixel, numberPixel);
        else
            for d = 1:DGuid
                M{d} = speye(numberPixel, numberPixel);
                MT{d} = speye(numberPixel, numberPixel);
            end
        end
    end
    
    f_old_vec = f_new_vec;
    
    for innerIter = 1:params.maxInnerIter
        
        
        %update the intermediate image f
        [f_new_vec,sigma] = update_fMC( f_new_vec,g,params.alpha,params.beta,...
            v_vec,bv_vec, bu_vec,u_vec,M,MT,grad,L,params.tol,...
            params.maxCgIter,Down,DownT,Dim,DGuid,weightColor,params.factor,params);
        if nargout > 1
            E = [E,calculateEnergyFunctionMC(f_old_vec, g, grad, sigma,...
                params.mu, params.lambda,params.patchsize, params.quantile,...
                params.guide,params.weighting,params.dev,params.cases,Down,...
                params.factor,weightColor, M)];
        end
        v_old = v_vec;
        % part for the TV prior
        gradient = cell(1,params.directions);
        if params.mu > 0
            % Gradient
            for k = 1:Dim
                start = (k-1)*numberPixel+1;
                ende = k*numberPixel;
                %Anisotropic TV denoising
                for pos = 1:params.directions
                    gradient{pos} = getGradient(f_new_vec(start:ende),grad,pos);
                    v_vec{pos}(start:ende) = update_vAnisotropic(...
                        gradient{pos},bv_vec{pos}(start:ende),params.beta,params.mu);
                    %update bregman variable bv of the TV
                    bv_vec{pos}(start:ende) = bv_vec{pos}(start:ende) + ...
                        (gradient{pos} - v_vec{pos}(start:ende));
                end
            end
        end
        u_old = u_vec;
        %part for the (A)quasi prior
        if params.lambda > 0
            for k = 1:Dim
                start = (k-1)*numberPixel+1;
                ende = k*numberPixel;
                if multiGuidance
                    u_vec(start:ende) = update_u( ...
                        weightColor(k)*f_new_vec(start:ende),M{k},...
                        bu_vec(start:ende),params.lambda,params.alpha);
                    bu_vec(start:ende) = update_bu(...
                        bu_vec(start:ende),weightColor(k)*f_new_vec(start:ende),...
                        u_vec(start:ende),M{k});
                else
                    u_vec(start:ende) = update_u( ...
                        weightColor(k)*f_new_vec(start:ende),M{1},...
                        bu_vec(start:ende),params.lambda,params.alpha);
                    bu_vec(start:ende) = update_bu( ...
                        bu_vec(start:ende),weightColor(k)*f_new_vec(start:ende),...
                        u_vec(start:ende),M{1});
                end
                %update the auxiliary variable u and the bregman variable bu
                
            end
        end
        

    end
    
    
    
    %compute the difference between the intermediate image and the
    %image from the previous iteration step
    %return is difference is small
    if params.stopEarly  && ((norm(f_old_vec - f_new_vec) / norm(f_old_vec)) < params.tol)
        f = vectorToMatrix(f_new_vec,R,C,Dim);
        
        return;
    end
    
end

f = vectorToMatrix(f_new_vec,R,C,Dim);