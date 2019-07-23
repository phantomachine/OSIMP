%Figure A1: impuse responses under cp vs. np
clear all;

pflag = 1; % = 1; if print eps files

n1 = 5;
nx = 5;
ymax = 3.0;

k = 5;

filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f')];

% home markup shock
uvec = zeros(21,1);
uvec(1) = 1.0;

eval(['load ./cp/mat/ifcp_', filename, '_rho1.0.mat']);
x1c_1 = y1vec - uvec;
x2c_1 = y2vec;
p1c_1 = p1vec;
p2c_1 = p2vec;
v1c_1 = v1vec;
v2c_1 = v2vec;
ttc_1 = ttvec;

eval(['load ./cp/mat/ifcp_', filename, '_rho2.0.mat']);
x1c_2 = y1vec - uvec;
x2c_2 = y2vec;
p1c_2 = p1vec;
p2c_2 = p2vec;
v1c_2 = v1vec;
v2c_2 = v2vec;
ttc_2 = ttvec;

eval(['load ./cp/mat/ifcp_', filename, '_rho0.5.mat']);
x1c_3 = y1vec - uvec;
x2c_3 = y2vec;
p1c_3 = p1vec;
p2c_3 = p2vec;
v1c_3 = v1vec;
v2c_3 = v2vec;
ttc_3 = ttvec;

eval(['load ./ncp/mat/ifnc_', filename, '_rho1.0.mat']);
x1n_1 = y1vec - uvec;
x2n_1 = y2vec;
p1n_1 = p1vec;
p2n_1 = p2vec;
v1n_1 = v1vec;
v2n_1 = v2vec;
ttn_1 = ttvec;

eval(['load ./ncp/mat/ifnc_', filename, '_rho2.0.mat']);
x1n_2 = y1vec - uvec;
x2n_2 = y2vec;
p1n_2 = p1vec;
p2n_2 = p2vec;
v1n_2 = v1vec;
v2n_2 = v2vec;
ttn_2 = ttvec;

eval(['load ./ncp/mat/ifnc_', filename, '_rho0.5.mat']);
x1n_3 = y1vec - uvec;
x2n_3 = y2vec;
p1n_3 = p1vec;
p2n_3 = p2vec;
v1n_3 = v1vec;
v2n_3 = v2vec;
ttn_3 = ttvec;


T = 4;

figure;
subplot(5,3,1)
plot([0:T-1],x1c_3(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],x1n_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-2.5 0.5]);
ylim([-1.25 0.05]);
ylabel('H gap');
title('\rho<1');

subplot(5,3,2)
plot([0:T-1],x1c_1(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],x1n_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-2.5 0.5]);
ylim([-1.25 0.05]);
title('\rho=1');

subplot(5,3,3)
plot([0:T-1],x1c_2(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],x1n_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-2.5 0.5]);
ylim([-1.25 0.05]);
title('\rho>1');
legend('CP','NP','Location','SouthEast');

subplot(5,3,4)
plot([0:T-1],p1c_3(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],p1n_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.05 0.05]);
ylim([-0.02 0.035]);
ylabel('H inf');

subplot(5,3,5)
plot([0:T-1],p1c_1(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],p1n_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.05 0.05]);
ylim([-0.02 0.035]);

subplot(5,3,6)
plot([0:T-1],p1c_2(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],p1n_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.05 0.05]);
ylim([-0.02 0.035]);

subplot(5,3,7)
plot([0:T-1],x2c_3(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],x2n_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.3 0.3]);
ylim([-0.15 0.075]);
yticks([-0.15 0]);
ylabel('F gap');

subplot(5,3,8)
plot([0:T-1],x2c_1(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],x2n_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.3 0.3]);
ylim([-0.15 0.075]);
yticks([-0.15 0]);

subplot(5,3,9)
plot([0:T-1],x2c_2(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],x2n_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.3 0.3]);
ylim([-0.15 0.075]);
yticks([-0.15 0]);

subplot(5,3,10)
plot([0:T-1],p2c_3(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],p2n_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.015 0.015]);
ylim([-0.01 0.01]);
yticks([-0.01 0 0.01]);
ylabel('F inf');

subplot(5,3,11)
plot([0:T-1],p2c_1(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],p2n_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.015 0.015]);
ylim([-0.01 0.01]);
yticks([-0.01 0 0.01]);

subplot(5,3,12)
plot([0:T-1],p2c_2(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],p2n_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.015 0.015]);
ylim([-0.01 0.01]);
yticks([-0.01 0 0.01]);

subplot(5,3,13)
plot([0:T-1],ttc_3(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],ttn_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
ylim([-0.25 0.625]);
yticks([0 0.5]);
ylabel('ToT');

subplot(5,3,14)
plot([0:T-1],ttc_1(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],ttn_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
ylim([-0.25 0.625]);
yticks([0 0.5]);

subplot(5,3,15)
plot([0:T-1],ttc_2(1:T),'g--','LineWidth',1.0);
hold on;
plot([0:T-1],ttn_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
ylim([-0.25 0.625]);
yticks([0 0.5]);

if(pflag); print -depsc2 figureA1_impulse_responses_cpandnp.eps; end;
