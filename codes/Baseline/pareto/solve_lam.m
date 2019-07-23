clear all;

% parameters are from Benigno and Benigno (2006) p.485
bet = 0.99;
eta = 0.47;
alp1 = 0.66; % for foreign country in the original calibration
alp2 = alp1;
rho = 1.0;
sig = 10.0;  
gam = 1/2;
kap1 = (1-alp1)*(1-alp1*bet)/alp1/(1+sig*eta);
kap2 = (1-alp2)*(1-alp2*bet)/alp2/(1+sig*eta);
sig1 = 1.0;
sig2 = 1/5;

ymax = 3.0;
n1 = 5;

knotsy = linspace(-ymax,ymax,n1)';
ry = n1-2;
invTy = inv(spbas(ry,knotsy));

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


k = 5;
nx = 5;
midy = ceil(n1/2);
midx = ceil(nx/2);

rvec = linspace(0.5,3.0,6);

for j = 1:size(rvec,2)

    rho = rvec(j)
%    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];

    % load values under scp (\nu=1/2)
    eval(['load ../scp/mat/pfsc_', filename, '.mat  v1mat v2mat']);  
%    eval(['load ../scp/mat/pfsc_k5_ny5nx5ymax3.0_rho', num2str(rho,'%1.1f'), '.mat v1mat v2mat']);
    v1s = v1mat(5,5,midy,midy,midx);
    v2s = v2mat(5,5,midy,midy,midx);

    fprintf('solving for the pareto weight minimizing the distance between vp and vs\n');
    [lam rc] = csolve('pareto',0.5,[],1e-6,100,Gu,Pu,invTy,bet,eta,rho,sig,kap1,kap2,gam,knotsy,v1s,v2s)
    fprintf('\n');

    [y1mat y2mat p1mat p2mat v1mat v2mat] = calcpar(Gu,Pu,invTy,invTy,bet,eta,rho,sig,kap1,kap2,gam,lam,knotsy,knotsy);
    plotif

    eval(['save ./mat/ifpar_', filename, '.mat y1vec y2vec p1vec p2vec v1vec v2vec ttvec lam;']);    
%    eval(['save ./mat/ifpar_k5_ny5nx5ymax3.0_rho', num2str(rho,'%1.1f'), '.mat y1vec y2vec p1vec p2vec v1vec v2vec ttvec lam']);
    
end