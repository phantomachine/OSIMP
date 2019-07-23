%Figure 7: impulse responses under cp, np, scp, and pareto
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

eval(['load ./pareto/mat/ifpar_', filename_scp, '_rho1.0.mat']);
x1p_1 = y1vec - uvec;
x2p_1 = y2vec;
p1p_1 = p1vec;
p2p_1 = p2vec;
v1p_1 = v1vec;
v2p_1 = v2vec;
ttp_1 = ttvec;
lam_1 = lam;

eval(['load ./pareto/mat/ifpar_', filename_scp, '_rho2.0.mat']);
x1p_2 = y1vec - uvec;
x2p_2 = y2vec;
p1p_2 = p1vec;
p2p_2 = p2vec;
v1p_2 = v1vec;
v2p_2 = v2vec;
ttp_2 = ttvec;
lam_2 = lam;

eval(['load ./pareto/mat/ifpar_', filename_scp, '_rho0.5.mat']);
x1p_3 = y1vec - uvec;
x2p_3 = y2vec;
p1p_3 = p1vec;
p2p_3 = p2vec;
v1p_3 = v1vec;
v2p_3 = v2vec;
ttp_3 = ttvec;
lam_3 = lam;


eval(['load ./scp/mat/ifsc_', filename_scp, '_rho1.0.mat']);
x1s_1 = y1vec - uvec;
x2s_1 = y2vec;
p1s_1 = p1vec;
p2s_1 = p2vec;
v1s_1 = v1vec;
v2s_1 = v2vec;
tts_1 = ttvec;
zs_1 = zvec;
xs_1 = xvec;

eval(['load ./scp/mat/ifsc_', filename_scp, '_rho2.0.mat']);
x1s_2 = y1vec - uvec;
x2s_2 = y2vec;
p1s_2 = p1vec;
p2s_2 = p2vec;
v1s_2 = v1vec;
v2s_2 = v2vec;
tts_2 = ttvec;
zs_2 = zvec;
xs_2 = xvec;

eval(['load ./scp/mat/ifsc_', filename_scp, '_rho0.5.mat']);
x1s_3 = y1vec - uvec;
x2s_3 = y2vec;
p1s_3 = p1vec;
p2s_3 = p2vec;
v1s_3 = v1vec;
v2s_3 = v2vec;
tts_3 = ttvec;
zs_3 = zvec;
xs_3 = xvec;

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


T = 7;

figure;
subplot(3,1,1)
plot([0:T-1],tts_3(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],ttp_3(1:T),'m:','LineWidth',2.0);
plot([0:T-1],ttc_3(1:T),'g--','LineWidth',1.0);
plot([0:T-1],ttn_3(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
ylim([-0.25 0.625]);
yticks([-0.25 0 0.5]);
title('Terms of Trade');
ylabel('\rho<1');
legend('SCP','Pareto','CP','NP','Location','NorthEast');

subplot(3,1,2)
plot([0:T-1],tts_1(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],ttp_1(1:T),'m:','LineWidth',2.0);
plot([0:T-1],ttc_1(1:T),'g--','LineWidth',1.0);
plot([0:T-1],ttn_1(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
ylim([-0.25 0.625]);
yticks([-0.25 0 0.5]);
ylabel('\rho=1');

subplot(3,1,3)
plot([0:T-1],tts_2(1:T),'b-','LineWidth',2.0);
hold on;
plot([0:T-1],ttp_2(1:T),'m:','LineWidth',2.0);
plot([0:T-1],ttc_2(1:T),'g--','LineWidth',1.0);
plot([0:T-1],ttn_2(1:T),'r-.','LineWidth',1.0);
xlim([0 T-1]);
ylim([-0.25 0.625]);
yticks([-0.25 0 0.5]);
ylabel('\rho>1');
%xlabel('Time');

if(pflag); print -depsc2 figure7_impulse_responses_pareto.eps; end;