% Figure 1: welfare comparison under cp vs. np
clear all;

pflag = 1; % = 1; if print eps files

n1 = 5;
ymax = 3.0;

nr = 46;
rvec = linspace(0.5,5.0,nr);
r0mat = zeros(nr,2,4);

for ir = 1:nr
    
    rho = rvec(ir);

    k = 5;
%     g = 4;
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
%    filename = ['g', num2str(g,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
    eval(['load ./cp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vc1 = ev1(3);
    vc2 = ev1(4);
    eval(['load ./ncp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vn1 = ev1(3);
    vn2 = ev1(4);
    r0mat(ir,1,1) = -max(vc1./vn1-1,0);
    r0mat(ir,2,1) = -max(vc2./vn2-1,0);
    
%     k = 4;
%     filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
%     rho = rvec(i);
%     eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec'])
%     x0vec(i,2) = xvec(20);
%     z0vec(i,2) = zvec(2);
%     
%     k = 3;
%     filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
%     rho = rvec(i);
%     eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec'])
%     x0vec(i,3) = xvec(20);
%     z0vec(i,3) = zvec(2);
     
    k = 2;
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    rho = rvec(ir);
    eval(['load ./cp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vc1 = ev1(3);
    vc2 = ev1(4);
    eval(['load ./ncp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vn1 = ev1(3);
    vn2 = ev1(4);
    r0mat(ir,1,2) = -max(vc1./vn1-1,0);
    r0mat(ir,2,2) = -max(vc2./vn2-1,0);
    
end


figure;
subplot(211);
plot(rvec,r0mat(:,1,1),'k-','LineWidth',3.0,'Color',[.5 .5 .5]);
hold on;
plot(rvec,r0mat(:,1,2),'k:','LineWidth',3.0,'Color',[.5 .5 .5]);
legend('\sigma/\sigma^*=5','\sigma/\sigma^*=2','Location','SouthEast')
title('-max(V^c/V^n-1,0)')
ylabel('Home','FontSize',16);
subplot(212);
plot(rvec,r0mat(:,2,1),'k-','LineWidth',3.0,'Color',[.5 .5 .5]);
hold on;
plot(rvec,r0mat(:,2,2),'k:','LineWidth',3.0,'Color',[.5 .5 .5]);
title('-max(V^{*c}/V^{*n}-1,0)')
xlabel('\rho','FontSize',16);
ylabel('Foreign','FontSize',16);

if(pflag); print -depsc2 figure1_welfare_comparison.eps; end;
