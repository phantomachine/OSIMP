function [ev1_decomp ev2_decomp] = decompev(rvec,sdmat,bet,eta,sig,kap1,kap2,gam)

nr = size(rvec,1);

for ir = 1:nr

    rho = rvec(ir);
    % for parameters: eta, rho, gam, sig, kap1, kap2,
%    eval(['load ./cp/mat/pfcp_', filename, '_rho', num2str(rho,'%1.1f'), '.mat'])
    
    sd = reshape(sdmat(:,ir),[7 1]);
    sd1 = zeros(5,1);
    sd2 = zeros(5,1);
    sd1(1) = gam*(eta+rho)*sd(1)^2;
    sd1(2) = (1-gam)*(eta+rho)*sd(2)^2;
    sd1(3) = gam*(1-gam)*(1-rho)*sd(3)^2;
    sd1(4) = gam*sig/kap1*sd(4)^2;
    sd1(5) = (1-gam)*sig/kap2*sd(5)^2;
    
    sd2(1) = (1-gam)*(eta+rho)*sd(6)^2;
    sd2(2) = gam*(eta+rho)*sd(7)^2;
    sd2(3:5) = sd1(3:5);

    ev1_decomp(6,ir) = -sum(sd1)/(1-bet);
    ev1_decomp(1:5,ir) = -sd1(1:5)/(1-bet);
    ev2_decomp(6,ir) = -sum(sd2)/(1-bet);
    ev2_decomp(1:5,ir) = -sd2(1:5)/(1-bet);

end