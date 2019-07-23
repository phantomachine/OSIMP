% Figure 2 and 5: welfare decomposition
clear all;

pflag = 1; % = 1; if print eps files

n1 = 5;
nx = 5;
ymax = 3.0;

nr = 46;
rvec = linspace(.5,5.0,nr)';
sdmat = zeros(3,7,nr);
evmat = zeros(3,4,nr);

k = 5;

for ir = 1:nr

    rho = rvec(ir);
    filename_scp = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];

    eval(['load ./cp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    sdmat(1,:,ir) = sd1;
    evmat(1,:,ir) = ev1;
    eval(['load ./ncp/sim/simresult_', filename, '.mat sd1 ev1 er1']);
    sdmat(2,:,ir) = sd1;
    evmat(2,:,ir) = ev1;
    eval(['load ./scp/sim/simresult_', filename_scp, '.mat sd1 ev1 er1 xx zz']);
    sdmat(3,:,ir) = sd1;
    evmat(3,:,ir) = ev1;
    
end

% load parameters except for rho
eval(['load ./cp/mat/pfcp_', filename, '.mat bet eta sig kap1 kap2 gam'])

% cooperation
[ev1_decomp ev2_decomp] = decompev(rvec,reshape(sdmat(1,:,:),[7 nr]),bet,eta,sig,kap1,kap2,gam);

figure;
subplot(211);
plot(rvec,reshape(evmat(1,1,:),[nr 1]),'g-','LineWidth',2.0);
bar2stack(rvec,ev1_decomp(1:5,:)');
hold on;
plot(rvec,reshape(evmat(1,1,:),[nr 1]),'g-','LineWidth',2.0);
xlim([rvec(1) rvec(end)]);
ylim([-40 10]);
legend('CP','H gap','F gap*','ToT','H inf','F inf','Location','SouthEast');
ylabel('Home','FontSize',16);
subplot(212);
plot(rvec,reshape(evmat(1,2,:),[nr 1]),'g-','LineWidth',2.0);
bar2stack(rvec,ev2_decomp(1:5,:)');
hold on;
plot(rvec,reshape(evmat(1,2,:),[nr 1]),'g-','LineWidth',2.0);
xlim([rvec(1) rvec(end)]);
ylim([-80 10]);
ylabel('Foreign','FontSize',16);
xlabel('\rho','FontSize',16);
legend('CP','F gap','H gap*','ToT','H inf','F inf','Location','SouthEast');

if(pflag); print -depsc2 figure2_welfare_decomposition_cp.eps; end;


% non-cooperation
[ev1_decomp ev2_decomp] = decompev(rvec,reshape(sdmat(2,:,:),[7 nr]),bet,eta,sig,kap1,kap2,gam);

figure;
subplot(211);
plot(rvec,reshape(evmat(2,1,:),[nr 1]),'r-','LineWidth',2.0);
bar2stack(rvec,ev1_decomp(1:5,:)');
hold on;
plot(rvec,reshape(evmat(2,1,:),[nr 1]),'r-','LineWidth',2.0);
xlim([rvec(1) rvec(end)]);
ylim([-40 10]);
legend('NCP','H gap','F gap*','ToT','H inf','F inf','Location','SouthEast');
ylabel('Home','FontSize',16);
subplot(212);
plot(rvec,reshape(evmat(2,2,:),[nr 1]),'r-','LineWidth',2.0);
bar2stack(rvec,ev2_decomp(1:5,:)');
hold on;
plot(rvec,reshape(evmat(2,2,:),[nr 1]),'r-','LineWidth',2.0);
xlim([rvec(1) rvec(end)]);
ylim([-80 10]);
ylabel('Foreign','FontSize',16);
xlabel('\rho','FontSize',16);
legend('NCP','F gap','H gap*','ToT','H inf','F inf','Location','SouthEast');

if(pflag); print -depsc2 figure2_welfare_decomposition_np.eps; end;


% sustainable cooperation
[ev1_decomp ev2_decomp] = decompev(rvec,reshape(sdmat(3,:,:),[7 nr]),bet,eta,sig,kap1,kap2,gam);

figure;
subplot(211);
plot(rvec,reshape(evmat(3,1,:),[nr 1]),'b-','LineWidth',2.0);
bar2stack(rvec,ev1_decomp(1:5,:)');
hold on;
plot(rvec,reshape(evmat(3,1,:),[nr 1]),'b-','LineWidth',2.0);
xlim([rvec(1) rvec(end)]);
ylim([-40 10]);
legend('SCP','H gap','F gap*','ToT','H inf','F inf','Location','SouthEast');
ylabel('Home','FontSize',16);
subplot(212);
plot(rvec,reshape(evmat(3,2,:),[nr 1]),'b-','LineWidth',2.0);
bar2stack(rvec,ev2_decomp(1:5,:)');
hold on;
plot(rvec,reshape(evmat(3,2,:),[nr 1]),'b-','LineWidth',2.0);
xlim([rvec(1) rvec(end)]);
ylim([-80 10]);
ylabel('Foreign','FontSize',16);
xlabel('\rho','FontSize',16);
legend('SCP','F gap','H gap*','ToT','H inf','F inf','Location','SouthEast');

if(pflag); print -depsc2 figure5_welfare_decomposition_scp.eps; end;
