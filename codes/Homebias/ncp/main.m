%clear all;

% parameters are from Benigno and Benigno (2006) p.485
bet = 0.99;
eta = 0.47;
alp = 0.66; %// 0.75 for foreign country in the original calibration
rho = 1.0;
sig = 10.0; 
neu = 0.6;
k1 = (1-alp)*(1-alp*bet)/alp/(1+sig*eta);
k2 = (1-alp)*(1-alp*bet)/alp/(1+sig*eta);

sig1 = 1.0;
sig2 = 1.0;
ymax = 3.0;
nu = 9;
n1 = 5;
n2 = 5;
%nx = 7;

% [Gu Pu] = Rouwenhorst(0.0,1.0/sqrt(nu-1),nu);
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

knots1 = linspace(-ymax,ymax,n1)';
knots2 = linspace(-ymax,ymax,n2)';
knotsy = knots1;

tic;

% [Gu Pu] = Rouwenhorst(0.0,1.0/sqrt(nu-1),nu);
% knotsy = linspace(-5.0,5.0,ny)';
ry = n1-2;
invTy = inv(spbas(ry,knotsy));

[y1mat y2mat p1mat p2mat v1mat v2mat] = calcnc(Gu,Pu,invTy,invTy,knotsy,knotsy,bet,eta,rho,sig,k1,k2,neu);

eval(['save ./pfnc_k3_ny5nx5rho1.0_neu', num2str(neu,'%1.1f'), '.mat bet eta rho sig k1 k2 neu Gu Pu knots1 knots2 y1mat y2mat p1mat p2mat v1mat v2mat;']);    
