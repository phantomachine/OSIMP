% Figure 11: asymptotic pseudo weight, home bias (in Sec 4.3)
clear all;

pflag = 1; % = 1; if print eps files

n1 = 5;
nx = 5;
ymax = 3.0;

nr = 41;
rvec = linspace(0.8,1.2,nr);
rvec(40) = 1.2; % 1.19?
x0mat = zeros(nr,4);
z0mat = zeros(nr,4);

for ir = 1:nr

    neu = rvec(ir);

    k = 3;
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
    eval(['load ./scp/mat/ifsc_', filename, '_neu', num2str(neu,'%1.2f'), '.mat xvec zvec']);    
    x0mat(ir,1) = xvec(20);
    z0mat(ir,1) = zvec(2);
    
end


figure;
plot(rvec,x0mat(:,1),'b-','LineWidth',3.0);
hold on;
plot([rvec(1) rvec(end)],[.5 .5],'r-');
xlabel('\rho','FontSize',16);
ylabel('asymptotic pseudo weight');
ylim([0.5 0.6]);

if(pflag); print -depsc2 figure11_pseudo_weight_neu.eps; end;
