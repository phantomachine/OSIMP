% Figure 4: asymptotic pseudo weight
clear all;

pflag = 1; % = 1; if print eps files

n1 = 5;
nx = 5;
ymax = 3.0;

nr = 46;
rvec = linspace(0.5,5.0,nr);
x0mat = zeros(nr,4);
z0mat = zeros(nr,4);

for ir = 1:nr

    rho = rvec(ir);

    k = 5;
%     g = 4;
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
%    filename = ['g', num2str(g,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
    eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec']);    
    x0mat(ir,1) = xvec(20);
    z0mat(ir,1) = zvec(2);
    
%     k = 4;
%     filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
%     eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec'])
%     x0vec(i,2) = xvec(20);
%     z0vec(i,2) = zvec(2);
%     
%     k = 3;
%     filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
%     eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec'])
%     x0vec(i,3) = xvec(20);
%     z0vec(i,3) = zvec(2);
     
    k = 2;
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
    eval(['load ./scp/mat/ifsc_', filename, '_rho', num2str(rho,'%1.1f'), '.mat xvec zvec']);    
    x0mat(ir,4) = xvec(20);
    z0mat(ir,4) = zvec(2);
    
end


figure;
plot(rvec,x0mat(:,1),'b-','LineWidth',3.0);
hold on;
%plot(rvec,x0vec(:,2),'m-.','LineWidth',2.0);
% plot(rvec,x0vec(:,3),'c--','LineWidth',2.0);
plot(rvec,x0mat(:,4),'b:','LineWidth',3.0);
plot([rvec(1) rvec(end)],[.5 .5],'r-');
legend('\sigma_{\epsilon}/\sigma^*_{\epsilon}=5','\sigma_{\epsilon}/\sigma^*_{\epsilon}=2','Location','NorthEast')
xlabel('\rho','FontSize',16);
ylabel('asymptotic pseudo weight');
ylim([0.5 0.8]);

if(pflag); print -depsc2 figure4_pseudo_weight.eps; end;
