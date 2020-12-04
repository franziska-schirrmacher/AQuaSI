function [ f,E] = admmSC(g,params)
%store the size of the image


[params,f_new_vec,Down,DownT,E,grad,...
    v_vec,bv_vec,u_vec,bu_vec,L] = setParamsSC(g,params);
[R,C] = size(params.guide);
numberPixel = R*C;

for outerIter = 1:params.maxOuterIter
    
    if params.lambda > 0
        % Compute look up table for quantile regularization (M = Id - Q)
        
        switch params.weighting
            case 'weighted'
                %%Compute the weighted median of the neighbourhoof
                switch params.cases
                    
                    case 'intermediate'
                        [LookUpRow, LookUpColumn] = main_aquasi_mex(...
                            vectorToMatrix(f_new_vec,R,C),...
                            vectorToMatrix(f_new_vec,R,C),...
                            params.patchsize,params.dev,params.quantile);
%                         [LookUpRow, LookUpColumn] = compute_Q(...
%                             vectorToMatrix(f_new_vec,R,C), ...
%                             params.patchsize,params.quantile,...
%                             vectorToMatrix(f_new_vec,R,C),...
%                             params.weighting,params.dev, ...
%                             params.cases,weightColor);
                    case 'guidance'
                        [LookUpRow, LookUpColumn] = main_aquasi_mex(...
                            vectorToMatrix(f_new_vec,R,C),...
                            params.guide,params.patchsize,params.dev,...
                            params.quantile);
%                         [LookUpRow, LookUpColumn] = compute_Q(...
%                             vectorToMatrix(f_new_vec,R,C), ...
%                             params.patchsize,params.quantile,...
%                             params.guide,...
%                             params.weighting,params.dev, ...
%                             params.cases,weightColor);

                end
            case 'unweighted'
                %%Compute the median of the neighbourhood
                [LookUpRow, LookUpColumn] = main_quasi_mex(...
                    vectorToMatrix(f_new_vec,R,C),params.guide,...
                    params.patchsize,params.dev,params.quantile);
%              [LookUpRow, LookUpColumn] = compute_Q(...
%                   vectorToMatrix(f_new_vec,R,C), ...
%                   params.patchsize,params.quantile,...
%                   params.guide,...
%                   params.weighting,params.dev, ...
%                   params.cases,weightColor);
        end
        
        %[LookUpRow, LookUpColumn] = compute_Q( vectorToMatrix(f_new_vec,R,C,Dim),patchsize,quantile,vectorToMatrix(f_new_vec,R,C,Dim),weighting,dev,guided,weightColor);
        M = speye(numberPixel, numberPixel) - sparse(double(LookUpRow),double(LookUpColumn),1,numberPixel,numberPixel);
        MT = speye(numberPixel, numberPixel) - sparse(double(LookUpColumn),double(LookUpRow),1,numberPixel,numberPixel);
    else
        M = speye(numberPixel, numberPixel);
        MT = speye(numberPixel, numberPixel);
    end
    
    f_old_vec = f_new_vec;
    
    for innerIter = 1:params.maxInnerIter
        
        %update the intermediate image f
        [f_new_vec,sigma] = update_f( f_new_vec,g,params.alpha,params.beta,v_vec,bv_vec,...
            bu_vec,u_vec,M,MT,grad,L,params.tol,params.maxCgIter,Down,DownT,params.factor,params);
        
        
        v_old = v_vec;
        % part for the TV prior
        gradient = cell(1,params.directions);
        if params.mu > 0
            % Gradient
            %Anisotropic TV denoising
            for pos = 1:params.directions
                gradient{pos} = getGradient(f_new_vec,grad,pos);
                v_vec{pos} = update_vAnisotropic(gradient{pos},bv_vec{pos},params.beta,params.mu);
                %update bregman variable bv of the TV
                bv_vec{pos} = bv_vec{pos} + (gradient{pos} - v_vec{pos});
            end
            
        end
        u_old = u_vec;
        %part for the (A)quasi prior
        if params.lambda > 0
            %update the auxiliary variable u and the bregman variable bu
            u_vec = update_u( f_new_vec,M,bu_vec,params.lambda,params.alpha);
            bu_vec = update_bu( bu_vec,f_new_vec,u_vec,M);
        end
        
        if nargout > 1
            E = [E,calculateEnergyFunction(f_new_vec, g, grad, sigma, params.mu, params.lambda,...
                params.patchsize, params.quantile,params.guide,params.weighting,params.dev,params.cases,Down,params.factor)];
        end
        %% Adapt alpha and beta
        
    end
    
   
 
    %compute the difference between the intermediate image and the
    %image from the previous iteration step
    %return is difference is small
    if  params.stopEarly  && (norm(f_old_vec - f_new_vec) / norm(f_old_vec) < params.tol)
        f = vectorToMatrix(f_new_vec,R,C);
        
        return;
    end
    
end

f = vectorToMatrix(f_new_vec,R,C);