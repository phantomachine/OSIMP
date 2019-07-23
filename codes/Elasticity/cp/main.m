%clear all;

% parameters are from Benigno and Benigno (2006) p.485
bet = 0.99;
eta = 0.47;
alp1 = 0.75; %0.66; % for foreign country in the original calibration
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

[y1mat y2mat p1mat p2mat v1mat v2mat] = calccp(Gu,Pu,invTy,invTy,bet,eta,rho,sig,kap1,kap2,gam,knotsy,knotsy);

plotif;

figure;
subplot(321);
plot([1:T],y1vec(1:end));
xlim([1 T]);
subplot(322);
plot([1:T],y2vec(1:end));
xlim([1 T]);
subplot(323);
plot([1:T],p1vec(1:end));
xlim([1 T]);
%ylim([-0.01 0.02]);
subplot(324);
plot([1:T],p2vec(1:end));
xlim([1 T]);
subplot(325);
plot([1:T],ttvec(1:end));
xlim([1 T]);
