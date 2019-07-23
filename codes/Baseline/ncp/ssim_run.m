clear all;

k = 5;
n1 = 5;
ymax = 3.0;

rvec = linspace(.5,5.0,46);

for i = 1:size(rvec,2)

    rho = rvec(i);
    filename = ['k', num2str(k,'%1d'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
    disp(filename);
    [sd1 ev1 er1] = stochsim(filename);
    eval(['save ./sim/simresult_', filename, '.mat sd1 ev1 er1']);
%    sdmat(i,:) = sd1;
    
end