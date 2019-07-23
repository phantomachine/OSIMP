% Figure 6: Pareto frontier
clear all;

pflag = 1; % = 1; if print eps files

k = 5;
n1 = 5;
nx = 5;
ymax = 3.0;

figure;
rvec = linspace(0.5,1.5,3);

for i = 1:size(rvec,2)

    rho = rvec(i);
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    disp(filename);
    eval(['load ./pareto/mat/pareto_', filename '.mat']);

    subplot(1,3,i);
    plot(1-v1s/v1n,1-v2s/v2n,'b+','MarkerSize',10.0);
    hold on;
    plot(1-v1s_ub/v1n,1-v2s_ub/v2n,'c^','MarkerSize',10.0);
    plot(1-v1pvec/v1n,1-v2pvec/v2n,'m-','LineWidth',2.0);
    plot(1-v1s/v1n,1-v2s/v2n,'b+','MarkerSize',10.0);
    plot(1-v1s_ub/v1n,1-v2s_ub/v2n,'c^','MarkerSize',10.0);
    plot(1-v1c/v1n,1-v2c/v2n,'g*','MarkerSize',10.0);
    plot(0,0,'ro','MarkerSize',10.0);
    plot([0 0],[-1d+4 1d+4],'r:');
    plot([-1d+4 1d+4],[0 0],'r:');
    title(sprintf('\\rho = %1.1f',rho));
    xlim([-3 0.5]);
    ylim([-.05 1.15]);
%    if (i==3); legend('\nu=1/2','\nu=u.b.','Location','SouthWest'); end;
    axis square;

end

if(pflag); print -depsc2 figure6_pareto_frontier1.eps; end;

figure;
rvec = linspace(2.0,3.0,3);

for i = 1:size(rvec,2)

    rho = rvec(i);
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    disp(filename);
    eval(['load ./pareto/mat/pareto_', filename '.mat']);

    subplot(1,3,i);
    plot(1-v1s/v1n,1-v2s/v2n,'b+','MarkerSize',10.0);
    hold on;
    plot(1-v1s_ub/v1n,1-v2s_ub/v2n,'c^','MarkerSize',10.0);
    plot(1-v1pvec/v1n,1-v2pvec/v2n,'m-','LineWidth',2.0);
    plot(1-v1s/v1n,1-v2s/v2n,'b+','MarkerSize',10.0);
    plot(1-v1s_ub/v1n,1-v2s_ub/v2n,'c^','MarkerSize',10.0);
    plot(1-v1c/v1n,1-v2c/v2n,'g*','MarkerSize',10.0);
    plot(0,0,'ro','MarkerSize',10.0);
    plot([0 0],[-1d+4 1d+4],'r:');
    plot([-1d+4 1d+4],[0 0],'r:');
    title(sprintf('\\rho = %1.1f',rho));
    xlim([-3 0.5]);
    ylim([-.05 1.15]);
    if (i==3); legend('\nu=1/2','\nu=u.b.','Location','SouthWest'); end;
    axis square;

end

if(pflag); print -depsc2 figure6_pareto_frontier2.eps; end;

% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'centimeters');
% set(gcf, 'PaperPosition', [0 0 16.5 5]);
% set(gcf, 'PaperPosition', [-1.7 -0.1 19.85 5.1]);
% set(gcf, 'PaperSize', [16.5 5]);
% %saveas(gcf,'pareto.pdf', 'pdf')