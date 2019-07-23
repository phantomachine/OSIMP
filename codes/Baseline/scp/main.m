%clear all;

damp = 1.0;

% parameters are from Benigno and Benigno (2006) p.485
bet = 0.99;
eta = 0.47;
alp1 = 0.66; %0.75; % for foreign country in the original calibration
alp2 = alp1;
rho = 1.0;
%rho = 3.0;
sig = 10.0;  
gam = 1/2;
kap1 = (1-alp1)*(1-alp1*bet)/alp1/(1+sig*eta);
kap2 = (1-alp2)*(1-alp2*bet)/alp2/(1+sig*eta);
sig1 = 1.0;
sig2 = 1/5;

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

tic;

[w1mat w2mat w11mat w12mat w21mat w22mat] = calcnc(Gu,Pu,invT1,invT2,bet,eta,rho,sig,kap1,kap2,gam,knots1,knots2);

p1mat0 = zeros(nu,nu,n1,n2,nx);
p2mat0 = zeros(nu,nu,n1,n2,nx);
v1mat0 = zeros(nu,nu,n1,n2,nx);
v2mat0 = zeros(nu,nu,n1,n2,nx);
x1mat0 = zeros(nu,nu,n1,n2,nx);
x2mat0 = zeros(nu,nu,n1,n2,nx);
y1mat = zeros(nu,nu,n1,n2,nx);
y2mat = zeros(nu,nu,n1,n2,nx);
zmat  = ones(nu,nu,n1,n2,nx);
xmat  = .5*ones(nu,nu,n1,n2,nx);

crit = 1e-5;
diff = 1e+4;
iter = 0;

while (diff>crit)
    
    [p1mat p2mat v1mat v2mat x1mat x2mat y1mat y2mat zmat xmat] = ...
    sus_iter(p1mat0,p2mat0,v1mat0,v2mat0,x1mat0,x2mat0,y1mat,y2mat,zmat,xmat,...
    w1mat,w2mat,w11mat,w12mat,w21mat,w22mat,Gu,Pu,invT1,invT2,invTx,bet,eta,rho,sig,kap1,kap2,gam,knots1,knots2,knotsx);
    
    diffg1 = max(max(max(max(max(abs(p1mat-p1mat0))))));
    diffg2 = max(max(max(max(max(abs(p2mat-p2mat0))))));
    diffg = max(diffg1,diffg2);
    diffv1 = max(max(max(max(max(abs(v1mat-v1mat0))))));
    diffv2 = max(max(max(max(max(abs(v2mat-v2mat0))))));
    diffv = max(diffv1,diffv2);
    diff = max(diffv,diffg);
    iter = iter + 1;
    disp(sprintf(' iteration %4d:    ||v1-v0|| = %10.5f,  ||g1-g0|| = %10.5f',iter,diffv,diffg));
    
    x1mat0 = damp*x1mat + (1-damp)*x1mat0;
    x2mat0 = damp*x2mat + (1-damp)*x2mat0;
    p1mat0 = damp*p1mat + (1-damp)*p1mat0;
    p2mat0 = damp*p2mat + (1-damp)*p2mat0;
    v1mat0 = damp*v1mat + (1-damp)*v1mat0;
    v2mat0 = damp*v2mat + (1-damp)*v2mat0;
    
end

toc;

%eval(['save ./mat/pfsc_k5_ny7nx7ymax3.0_rho', num2str(rho,'%1.1f'), '.mat bet eta rho sig kap1 kap2 gam Gu Pu knots1 knots2 knotsx x1mat x2mat y1mat y2mat p1mat p2mat v1mat v2mat w1mat w2mat zmat xmat;']);