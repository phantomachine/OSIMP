% Figure 10: asymptotic pseudo weight, non-unitary elasticity
clear all;

pflag = 1; % = 1; if print eps files

n1 = 5;
nx = 5;
ymax = 3.0;

nr = 41;
rvec = linspace(1.0,5.0,nr);
x0mat = zeros(nr,1);
z0mat = zeros(nr,1);

for ir = 1:nr

    iota = rvec(ir);

    k = 5;
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'nx', num2str(nx,'%1d'), 'ymax', num2str(ymax,'%1.1f')];
    eval(['load ./scp/mat/ifsc_', filename, '_iota', num2str(iota,'%1.1f'), '.mat xvec zvec']);    
    x0mat(ir,1) = xvec(20);
    z0mat(ir,1) = zvec(2);
    
end


figure;
plot(rvec,x0mat(:,1),'b-','LineWidth',3.0);
hold on;
plot([rvec(1) rvec(end)],[.5 .5],'r-');
xlabel('\iota','FontSize',16);
ylabel('asymptotic pseudo weight');
%ylim([0.5 0.8]);

if(pflag); print -depsc2 figure10_pseudo_weight_iota.eps; end;
