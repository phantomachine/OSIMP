clear all;

k = 5;
n1 = 5;
ymax = 3.0;

rvec = linspace(.5,5.0,46);
gvec = [3 3.25 3.5 3.75];

for j = 1:size(gvec,2)
%for k = 5:-1:2

%    gam = 1/gvec(j);

    for i = 1:size(rvec,2)

        rho = rvec(i);
        filename = ['g', num2str(gvec(j),'%1.3f'), '_ny', num2str(n1,'%1d'), 'ymax', num2str(ymax,'%1.1f'), '_rho', num2str(rho,'%1.1f')];
        disp(filename);
        [sd1 ev1 er1] = stochsim(filename);
        eval(['save ./sim/simresult_', filename, '.mat sd1 ev1 er1']);
    %    sdmat(i,:) = sd1;

    end
    
end