function T = spbas(r,knots)
% evaluate basis function of cubic spline on r+2 knots
% Feb 4 2010 Takeki Sunakawa

T = zeros(r,r);

for i = 1:r
    
    k = i+1; % index on knots
    
    dt0 = knots(k)-knots(k-1);
    dt1 = knots(k+1)-knots(k);
    
    if i==1
        
        T(i,i) = 2*(dt0+dt1)-(dt1-dt0^2/dt1);
        T(i,i+1) = dt0+dt0^2/dt1;    

    elseif i==r
        
        T(i,i-1) = dt1+dt1^2/dt0;
        T(i,i) = 2*(dt0+dt1)-(dt0-dt1^2/dt0);
        
    else
        
        T(i,i-1) = dt0;
        T(i,i) = 2*(dt0+dt1);
        T(i,i+1) = dt1; 
    
    end
        
end