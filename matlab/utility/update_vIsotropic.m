function [ dx, dy] = update_vIsotropic(  bx,by,f_new,lambda,R,C )

      
%UPDATE_B_ISOTROPICTV Summary of this function goes here
%   Detailed explanation goes here   
        fx = [f_new(1:R,2:C),zeros(R,1)];
        %f[r+1][c]
        fy = [f_new(2:R,1:C);zeros(1,C)];
        
        gx =  fx - f_new ;     
        gy =  fy - f_new;

        %compute s
        s = sqrt((gx+bx).^2 + (gy+by).^2);
        tmp = s - (1/lambda);
        %max operator 
        t = tmp < 0;       
        dx = (gx+bx)./s;
        dx = dx.*tmp;
        dy = (gy + by)./s;
        dy = dy.*tmp;
        dx(t) = 0;
        dy(t) = 0;
        
        
        %use shrink for the boundaries
        baseX = fx(R,:)- f_new(R,:) + bx(R,:);
        %baseX(R,C) = 0;        
        dx(R,:) = shrink(baseX,1/lambda);       
        
        baseY = fy(:,C) - f_new(:,C) + by(:,C);
        %baseY(R,C) = 0;
        dy(:,C) = shrink(baseY,1/lambda);
       
end

