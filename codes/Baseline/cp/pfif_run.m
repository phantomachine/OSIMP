clear all;

% parameters are from Benigno and Benigno (2006) p.485
bet = 0.99;
eta = 0.47;
alp1 = 0.66; % for foreign country in the original calibration
alp2 = alp1;
%rho = 1.0;
sig = 10.0;  
gam = 1/2;
kap1 = (1-alp1)*(1-alp1*bet)/alp1/(1+sig*eta);
kap2 = (1-alp2)*(1-alp2*bet)/alp2/(1+sig*eta);
sig1 = 1.0;
%sig2 = 1/5;

ymax = 3.0;
n1 = 5;
knotsy = linspace(-ymax,ymax,n1)';
ry = n1-2;
invTy = inv(spbas(ry,knotsy));

rvec = linspace(0.5,5.0,46);

for k = 5:-1:2

    sig2 = 1/k;
    
    Gu = zeros(9,2);
    Gu(1,:) = [-sig1 -sig2];
    Gu(2,:) = [-sig1 0];
    Gu(3,:) = [-sig1 sig2];
    Gu(4,:) = [0 -sig2];
    Gu(5,:) = [0 0];
    Gu(6,:) = [0 sig2];
    Gu(7,:) = [sig1 -sig2];
    Gu(8,:) = [sig1 0];
    Gu(9,:) = [sig1 sig2];
    Pu = [.25 .5 .25;
        .25 .5 .25;
        .25 .5 .25];
    Pu = kron(Pu,Pu);

    for i = 1:size(rvec,2)

        rho = rvec(i);
        filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
        disp(filename);
        [y1mat y2mat p1mat p2mat v1mat v2mat] = calccp(Gu,Pu,invTy,invTy,bet,eta,rho,sig,kap1,kap2,gam,knotsy,knotsy);
        eval(['save ./mat/pfcp_', filename, '.mat bet eta rho sig kap1 kap2 gam Gu Pu knotsy y1mat y2mat p1mat p2mat v1mat v2mat;']);    
%        eval(['load ./mat/pfcp_', filename, '.mat bet eta rho sig kap1 kap2 gam Gu Pu knotsy y1mat y2mat p1mat p2mat v1mat v2mat;']);    
        plotif;
        eval(['save ./mat/ifcp_', filename, '.mat y1vec y2vec p1vec p2vec v1vec v2vec ttvec;']);    

    end
    
end