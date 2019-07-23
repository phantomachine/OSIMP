% Figure 8: welfare comparison under cp vs. np, asymmetric country size (in
% Sec. 4.1)
clear all;

pflag = 1; % = 1; if print eps files

n1 = 5;
ymax = 3.0;

nr = 46;
rvec = linspace(0.5,5.0,nr);
r0mat = zeros(nr,2,4);

for ir = 1:nr
    
    rho = rvec(ir);

    g = 3.75;
    filename = ['g', num2str(g,'%1.3f'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    eval(['load ./cp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vc1 = ev1(3);
    vc2 = ev1(4);
    eval(['load ./ncp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vn1 = ev1(3);
    vn2 = ev1(4);
    r0mat(ir,1,1) = -max(vc1./vn1-1,0);
    r0mat(ir,2,1) = -max(vc2./vn2-1,0);
    
    g = 3.5;
    filename = ['g', num2str(g,'%1.3f'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    eval(['load ./cp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vc1 = ev1(3);
    vc2 = ev1(4);
    eval(['load ./ncp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vn1 = ev1(3);
    vn2 = ev1(4);
    r0mat(ir,1,2) = -max(vc1./vn1-1,0);
    r0mat(ir,2,2) = -max(vc2./vn2-1,0);
    
    g = 3.25;
    filename = ['g', num2str(g,'%1.3f'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    eval(['load ./cp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vc1 = ev1(3);
    vc2 = ev1(4);
    eval(['load ./ncp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vn1 = ev1(3);
    vn2 = ev1(4);
    r0mat(ir,1,3) = -max(vc1./vn1-1,0);
    r0mat(ir,2,3) = -max(vc2./vn2-1,0);
     
    g = 3;
    filename = ['g', num2str(g,'%1.3f'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    eval(['load ./cp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vc1 = ev1(3);
    vc2 = ev1(4);
    eval(['load ./ncp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    vn1 = ev1(3);
    vn2 = ev1(4);
    r0mat(ir,1,4) = -max(vc1./vn1-1,0);
    r0mat(ir,2,4) = -max(vc2./vn2-1,0);
    
end


figure;
subplot(211);
plot(rvec,r0mat(:,1,4),'k:','LineWidth',3.0,'Color',[.5 .5 .5]);
hold on;
plot(rvec,r0mat(:,1,3),'k--','LineWidth',3.0,'Color',[.5 .5 .5]);
plot(rvec,r0mat(:,1,2),'k-.','LineWidth',3.0,'Color',[.5 .5 .5]);
plot(rvec,r0mat(:,1,1),'k-','LineWidth',3.0,'Color',[.5 .5 .5]);
title('-max(V^c/V^n-1,0)');
xlim([rvec(1) rvec(end)]);
ylabel('Home','FontSize',16);
subplot(212);
plot(rvec,r0mat(:,2,4),'k:','LineWidth',3.0,'Color',[.5 .5 .5]);
hold on;
plot(rvec,r0mat(:,2,3),'k--','LineWidth',3.0,'Color',[.5 .5 .5]);
plot(rvec,r0mat(:,2,2),'k-.','LineWidth',3.0,'Color',[.5 .5 .5]);
plot(rvec,r0mat(:,2,1),'k-','LineWidth',3.0,'Color',[.5 .5 .5]);
legend('1/\gamma=3','1/\gamma=3.25','1/\gamma=3.5','1/\gamma=3.75','Location','NorthEast');
title('-max(V^{*c}/V^{*n}-1,0)');
xlim([rvec(1) rvec(end)]);
xlabel('\rho','FontSize',16);
ylabel('Foreign','FontSize',16);

if(pflag); print -depsc2 figure8_welfare_comparison_gam.eps; end;
