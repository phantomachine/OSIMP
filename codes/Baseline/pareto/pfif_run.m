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
lamvec = linspace(0,1,21);

for j = 1:size(rvec,2)

    rho = rvec(j);
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    filename_scp = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    disp(filename_scp);

    for i = 1:size(lamvec,2)

        lam = lamvec(i)
        eval(['load ./mat/pfpar_', filename, '_lam', num2str(lam,'%1.2f'), '.mat bet eta rho sig kap1 kap2 gam lam Gu Pu knotsy y1mat y2mat p1mat p2mat v1mat v2mat']);  
        [y1mat y2mat p1mat p2mat v1mat v2mat] = calcpar(Gu,Pu,invTy,invTy,bet,eta,rho,sig,kap1,kap2,gam,lam,knotsy,knotsy);
        eval(['save ./mat/pfpar_', filename, '_lam', num2str(lam,'%1.2f'), '.mat bet eta rho sig kap1 kap2 gam lam Gu Pu knotsy y1mat y2mat p1mat p2mat v1mat v2mat']);  
        v1pvec(i) = v1mat(5,5,midy,midy);
        v2pvec(i) = v2mat(5,5,midy,midy);

    end
    
    % load conditional values under cp, np, and scp
    eval(['load ../ncp/mat/pfnc_', filename, '.mat v1mat v2mat']);  
    v1n = v1mat(5,5,midy,midy);
    v2n = v2mat(5,5,midy,midy);
    eval(['load ../cp/mat/pfcp_', filename, '.mat v1mat v2mat']);  
    v1c = v1mat(5,midy,midy);
    v2c = v2mat(5,midy,midy);    
    eval(['load ../scp/mat/pfsc_', filename_scp, '.mat knotsx v1mat v2mat']);  
    % \nu = 1/2
    v1s = v1mat(5,5,midy,midy,midx);
    v2s = v2mat(5,5,midy,midy,midx);
    % \nu = u.b.
    eval(['load ../scp/sim/simresult_', filename_scp, '.mat xx']);  
    v1svec = reshape(v1mat(5,5,midy,midy,:),[nx 1]);
    v2svec = reshape(v2mat(5,5,midy,midy,:),[nx 1]);
    % linear interpolation???
    v1s_ub = interp1(knotsx,v1svec,xx,'linear');
    v2s_ub = interp1(knotsx,v2svec,xx,'linear');
        
    eval(['save ./mat/pareto_', filename_scp, '.mat v1pvec v2pvec v1c v2c v1s v2s v1s_ub v2s_ub v1n v2n']);

%     figure;
%     plot(v1pvec,v2pvec,'m-','LineWidth',2.0);
%     hold on;
%     plot(v1c,v2c,'g*','MarkerSize',10.0);
%     plot(v1n,v2n,'ro','MarkerSize',10.0);
%     plot([v1n v1n],[-1d+4 1d+4],'r:');
%     plot([-1d+4 1d+4],[v2n v2n],'r:');
%     plot(v1s,v2s,'b+','MarkerSize',10.0);
%     xmin = min(v1pvec);
%     xmax = max(v1pvec);
%     ymin = min(v2pvec);
%     ymax1 = max(v2pvec);
%     zmin = min([xmin ymin]);
%     zmax = max([xmin ymax1]);
%     xlim([zmin zmax]);
%     ylim([zmin zmax]);
%     axis square;
%     drawnow;
    
end