clear all;

% parameters are from Benigno and Benigno (2006) p.485
bet = 0.99;
eta = 0.47;
alp1 = 0.66; % for foreign country in the original calibration
alp2 = alp1;
rho = 1.0;
sig = 10.0;  
%gam = 1/2;
kap1 = (1-alp1)*(1-alp1*bet)/alp1/(1+sig*eta);
kap2 = (1-alp2)*(1-alp2*bet)/alp2/(1+sig*eta);
sig1 = 1.0;
sig2 = 1/3;

ymax = 3.0;
nu = 9;
n1 = 5;
n2 = n1;
nx = 5;
r1 = n1-2;
rx = nx-2;
knots1 = linspace(-ymax,ymax,n1)';
knots2 = knots1;
knotsx = linspace(0.0,1.0,nx)';
knotsx(1) = 0.0001;
knotsx(end) = 0.9999;

invT1 = inv(spbas(r1,knots1));
invT2 = invT1;
invTx = inv(spbas(rx,knotsx));

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

k = 3;
rvec = linspace(0.8,1.2,41);
%rvec = [0.8 1.0 1.2];

for i = 1:size(rvec,2)

    neu = rvec(i);
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_neu', num2str(neu,'%1.2f')];
    disp(filename);
    calcscp;
    eval(['save ./mat/pfsc_', filename, '.mat bet eta rho sig kap1 kap2 neu Gu Pu knots1 knots2 knotsx x1mat x2mat y1mat y2mat p1mat p2mat v1mat v2mat zmat xmat w1mat w2mat;']);    
    plotif;
    eval(['save ./mat/ifsc_', filename, '.mat y1vec y2vec p1vec p2vec v1vec v2vec zvec xvec ttvec;']);    
    
end