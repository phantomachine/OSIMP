% Figure 9: asymptotic pseudo weight, asymmetric country size (in
% Sec. 4.1)
clear all;

pflag = 1; % = 1; if print eps files

n1 = 5;
nx = 5;
ymax = 3.0;

nr = 36;
rvec = linspace(0.5,4.0,nr);
x0mat = zeros(nr,4);
z0mat = zeros(nr,4);

for ir = 1:nr

    rho = rvec(ir);

    g = 3.75;
    filename = ['g', num2str(g,'%1.3f'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
    eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec']);    
    x0mat(ir,1) = xvec(20);
    z0mat(ir,1) = zvec(2);
    
    g = 3.5;
    filename = ['g', num2str(g,'%1.3f'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
    eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec']);    
    x0mat(ir,2) = xvec(20);
    z0mat(ir,2) = zvec(2);

    g = 3.25;
    filename = ['g', num2str(g,'%1.3f'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
    eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec']);    
    x0mat(ir,3) = xvec(20);
    z0mat(ir,3) = zvec(2);
    
    g = 3;
    filename = ['g', num2str(g,'%1.3f'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
    eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec']);    
    x0mat(ir,4) = xvec(20);
    z0mat(ir,4) = zvec(2);
    
end


figure;
plot(rvec,x0mat(:,4),'b:','LineWidth',3.0);
hold on;
plot(rvec,x0mat(:,3),'b--','LineWidth',3.0);
plot(rvec,x0mat(:,2),'b-.','LineWidth',3.0);
plot(rvec,x0mat(:,1),'b-','LineWidth',3.0);
plot([rvec(1) rvec(end)],[.5 .5],'r-');
legend('1/\gamma=3','1/\gamma=3.25','1/\gamma=3.5','1/\gamma=3.75','Location','NorthEast');
xlabel('\rho','FontSize',16);
ylabel('asymptotic pseudo weight');
ylim([0.2 0.6]);

if(pflag); print -depsc2 figure9_pseudo_weight_gam.eps; end;
