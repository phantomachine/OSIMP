function f = pareto(lam,Gu,Pu,invTy,bet,eta,rho,sig,kap1,kap2,gam,knotsy,v1s,v2s)


[y1mat y2mat p1mat p2mat v1mat v2mat] = calcpar(Gu,Pu,invTy,invTy,bet,eta,rho,sig,kap1,kap2,gam,lam,knotsy,knotsy);
fprintf('.');

v1p = v1mat(5,5,3,3);
v2p = v2mat(5,5,3,3);

%f = abs(v1p-v1s) + abs(v2p-v2s);
f = abs(v2p-v2s);