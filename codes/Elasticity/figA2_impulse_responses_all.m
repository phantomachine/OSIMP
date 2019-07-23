%Figure A2: impulse responses, non-unitary elasticity
clear all;

pflag = 1; % = 1; if print eps files

n1 = 5;
nx = 5;
ymax = 3.0;

k = 5;

filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
filename_scp = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
disp(filename_scp);

% home markup shock
uvec = zeros(21,1);
uvec(1) = 1.0;

eval(['load ./scp/mat/ifsc_', filename_scp, '_iota2.0.mat']);
x1s_1 = y1vec - uvec;
x2s_1 = y2vec;
p1s_1 = p1vec;
p2s_1 = p2vec;
v1s_1 = v1vec;
v2s_1 = v2vec;
tts_1 = ttvec;
zs_1 = zvec;
xs_1 = xvec;

eval(['load ./scp/mat/ifsc_', filename_scp, '_iota5.0.mat']);
x1s_2 = y1vec - uvec;
x2s_2 = y2vec;
p1s_2 = p1vec;
p2s_2 = p2vec;
v1s_2 = v1vec;
v2s_2 = v2vec;
tts_2 = ttvec;
zs_2 = zvec;
xs_2 = xvec;

eval(['load ./scp/mat/ifsc_', filename_scp, '_iota1.0.mat']);
x1s_3 = y1vec - uvec;
x2s_3 = y2vec;
p1s_3 = p1vec;
p2s_3 = p2vec;
v1s_3 = v1vec;
v2s_3 = v2vec;
tts_3 = ttvec;
zs_3 = zvec;
xs_3 = xvec;

eval(['load ./cp/mat/ifcp_', filename, '_iota2.0.mat']);
x1c_1 = y1vec - uvec;
x2c_1 = y2vec;
p1c_1 = p1vec;
p2c_1 = p2vec;
v1c_1 = v1vec;
v2c_1 = v2vec;
ttc_1 = ttvec;

eval(['load ./cp/mat/ifcp_', filename, '_iota5.0.mat']);
x1c_2 = y1vec - uvec;
x2c_2 = y2vec;
p1c_2 = p1vec;
p2c_2 = p2vec;
v1c_2 = v1vec;
v2c_2 = v2vec;
ttc_2 = ttvec;

eval(['load ./cp/mat/ifcp_', filename, '_iota1.0.mat']);
x1c_3 = y1vec - uvec;
x2c_3 = y2vec;
p1c_3 = p1vec;
p2c_3 = p2vec;
v1c_3 = v1vec;
v2c_3 = v2vec;
ttc_3 = ttvec;

eval(['load ./ncp/mat/ifnc_', filename, '_iota2.0.mat']);
x1n_1 = y1vec - uvec;
x2n_1 = y2vec;
p1n_1 = p1vec;
p2n_1 = p2vec;
v1n_1 = v1vec;
v2n_1 = v2vec;
ttn_1 = ttvec;

eval(['load ./ncp/mat/ifnc_', filename, '_iota5.0.mat']);
x1n_2 = y1vec - uvec;
x2n_2 = y2vec;
p1n_2 = p1vec;
p2n_2 = p2vec;
v1n_2 = v1vec;
v2n_2 = v2vec;
ttn_2 = ttvec;

eval(['load ./ncp/mat/ifnc_', filename, '_iota1.0.mat']);
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
plot([0:T-1],x1s_3(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],x1c_3(1:T),'g--','LineWidth',1.0);
plot([0:T-1],x1n_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-2.5 0.5]);
ylim([-1.25 0.05]);
ylabel('H gap');
title('\iota=1.0');

subplot(5,3,2)
plot([0:T-1],x1s_1(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],x1c_1(1:T),'g--','LineWidth',1.0);
plot([0:T-1],x1n_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-2.5 0.5]);
ylim([-1.25 0.05]);
title('\iota=2.0');

subplot(5,3,3)
plot([0:T-1],x1s_2(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],x1c_2(1:T),'g--','LineWidth',1.0);
plot([0:T-1],x1n_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-2.5 0.5]);
ylim([-1.25 0.05]);
title('\iota=5.0');
legend('SCP','CP','NP','Location','SouthEast');

subplot(5,3,4)
plot([0:T-1],p1s_3(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],p1c_3(1:T),'g--','LineWidth',1.0);
plot([0:T-1],p1n_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.05 0.05]);
ylim([-0.02 0.035]);
ylabel('H inf');

subplot(5,3,5)
plot([0:T-1],p1s_1(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],p1c_1(1:T),'g--','LineWidth',1.0);
plot([0:T-1],p1n_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.05 0.05]);
ylim([-0.02 0.035]);

subplot(5,3,6)
plot([0:T-1],p1s_2(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],p1c_2(1:T),'g--','LineWidth',1.0);
plot([0:T-1],p1n_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.05 0.05]);
ylim([-0.02 0.035]);

subplot(5,3,7)
plot([0:T-1],x2s_3(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],x2c_3(1:T),'g--','LineWidth',1.0);
plot([0:T-1],x2n_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.3 0.3]);
ylim([-0.15 0.3]);
yticks([-0.15 0]);
ylabel('F gap');

subplot(5,3,8)
plot([0:T-1],x2s_1(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],x2c_1(1:T),'g--','LineWidth',1.0);
plot([0:T-1],x2n_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.3 0.3]);
ylim([-0.15 0.3]);
yticks([-0.15 0]);

subplot(5,3,9)
plot([0:T-1],x2s_2(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],x2c_2(1:T),'g--','LineWidth',1.0);
plot([0:T-1],x2n_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.3 0.3]);
ylim([-0.15 0.3]);
yticks([-0.15 0]);

subplot(5,3,10)
plot([0:T-1],p2s_3(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],p2c_3(1:T),'g--','LineWidth',1.0);
plot([0:T-1],p2n_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.015 0.015]);
ylim([-0.01 0.01]);
yticks([-0.01 0 0.01]);
ylabel('F inf');

subplot(5,3,11)
plot([0:T-1],p2s_1(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],p2c_1(1:T),'g--','LineWidth',1.0);
plot([0:T-1],p2n_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.015 0.015]);
ylim([-0.01 0.01]);
yticks([-0.01 0 0.01]);

subplot(5,3,12)
plot([0:T-1],p2s_2(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],p2c_2(1:T),'g--','LineWidth',1.0);
plot([0:T-1],p2n_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
%ylim([-0.015 0.015]);
ylim([-0.01 0.01]);
yticks([-0.01 0 0.01]);

subplot(5,3,13)
plot([0:T-1],tts_3(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],ttc_3(1:T),'g--','LineWidth',1.0);
plot([0:T-1],ttn_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
ylim([-0.25 0.625]);
yticks([0 0.5]);
ylabel('ToT');

subplot(5,3,14)
plot([0:T-1],tts_1(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],ttc_1(1:T),'g--','LineWidth',1.0);
plot([0:T-1],ttn_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
ylim([-0.25 0.625]);
yticks([0 0.5]);

subplot(5,3,15)
plot([0:T-1],tts_2(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],ttc_2(1:T),'g--','LineWidth',1.0);
plot([0:T-1],ttn_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
ylim([-0.25 0.625]);
yticks([0 0.5]);

if(pflag); print -depsc2 figureA2_impulse_responses_iota.eps; end;
