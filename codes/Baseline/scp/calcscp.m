%clear all;

tic;
damp = 1.0;

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
    if(mod(iter,100)==0); disp(sprintf(' iteration %4d:    ||v1-v0|| = %10.5f,  ||g1-g0|| = %10.5f',iter,diffv,diffg)); end;
    
    x1mat0 = damp*x1mat + (1-damp)*x1mat0;
    x2mat0 = damp*x2mat + (1-damp)*x2mat0;
    p1mat0 = damp*p1mat + (1-damp)*p1mat0;
    p2mat0 = damp*p2mat + (1-damp)*p2mat0;
    v1mat0 = damp*v1mat + (1-damp)*v1mat0;
    v2mat0 = damp*v2mat + (1-damp)*v2mat0;

end

toc;