function [sd1 ev1 er1 ivec1 xvec1 zvec1] = stochsim(filename)

eval(['load ./mat/pfsc_', filename, '.mat'])

knotsy = knots1;
ny = length(knotsy);
nx = length(knotsx);
nu = length(Gu);
midy = ceil(ny/2);
midx = ceil(nx/2);

N = 1;
T = 2000;
drop = 1800; % last 200 periods are used
er_mean = zeros(4,N);
er_max = zeros(4,N);
sd = zeros(7,N);
ivec1 = zeros(T,N);
xvec1 = zeros(T,N);
zvec1 = zeros(T,N);

rng(0);

for i = 1:N

    % generate exogenous shock sequence
    drop = 0;
    ivec = zeros(T+1,1);
    ivec(1) = 5;
    cumP = cumsum(Pu')';

    for t = 1:T
        cumPi = cumP(ivec(t),:);
        ivec(t+1) = sum(rand-cumPi >= 0);
        ivec(t+1) = min(ivec(t+1)+1,nu);
    end

    y1vec = zeros(T+1,1);
    y2vec = zeros(T+1,1);
    p1vec = zeros(T+1,1);
    p2vec = zeros(T+1,1);
    v1vec = zeros(T+1,1);
    v2vec = zeros(T+1,1);
    w1vec = zeros(T+1,1);
    w2vec = zeros(T+1,1);
    xvec = gam*ones(T+1,1);
    zvec = ones(T+1,1);
    u1vec = zeros(T+1,1);
    u2vec = zeros(T+1,1);
    e1vec = zeros(T+1,1);
    e2vec = zeros(T+1,1);
    e3vec = zeros(T+1,1);
    e4vec = zeros(T+1,1);

    for t = 1:T

        ku = ivec(t);
        ju = ivec(t+1);

        y1vec(t+1) = intf3(knotsy,knotsy,knotsx,reshape(y1mat(ku,ju,:,:,:),[ny ny nx]),y1vec(t),y2vec(t),xvec(t));
        y2vec(t+1) = intf3(knotsy,knotsy,knotsx,reshape(y2mat(ku,ju,:,:,:),[ny ny nx]),y1vec(t),y2vec(t),xvec(t));
        p1vec(t+1) = intf3(knotsy,knotsy,knotsx,reshape(p1mat(ku,ju,:,:,:),[ny ny nx]),y1vec(t),y2vec(t),xvec(t));
        p2vec(t+1) = intf3(knotsy,knotsy,knotsx,reshape(p2mat(ku,ju,:,:,:),[ny ny nx]),y1vec(t),y2vec(t),xvec(t));    
        v1vec(t+1) = intf3(knotsy,knotsy,knotsx,reshape(v1mat(ku,ju,:,:,:),[ny ny nx]),y1vec(t),y2vec(t),xvec(t));
        v2vec(t+1) = intf3(knotsy,knotsy,knotsx,reshape(v2mat(ku,ju,:,:,:),[ny ny nx]),y1vec(t),y2vec(t),xvec(t));    
        w1vec(t+1) = intf2(knotsy,knotsy,reshape(w1mat(ku,ju,:,:),[ny ny]),y1vec(t),y2vec(t));
        w2vec(t+1) = intf2(knotsy,knotsy,reshape(w2mat(ku,ju,:,:),[ny ny]),y1vec(t),y2vec(t));    

        if ((v1vec(t+1)<=w1vec(t+1))||(v2vec(t+1)<=w2vec(t+1))) % approx error?
            zvec(t+1) = intf3(knotsy,knotsy,knotsx,reshape(zmat(ku,ju,:,:,:),[ny ny nx]),y1vec(t),y2vec(t),xvec(t));    
            xvec(t+1) = intf3(knotsy,knotsy,knotsx,reshape(xmat(ku,ju,:,:,:),[ny ny nx]),y1vec(t),y2vec(t),xvec(t));
        else
            zvec(t+1) = 1;
            xvec(t+1) = xvec(t);
        end

        u1vec(t+1) = Gu(ivec(t+1),1);
        u2vec(t+1) = Gu(ivec(t+1),2);

        y1 = y1vec(t+1);
        y2 = y2vec(t+1);
        x = xvec(t+1);
        p1 = p1vec(t+1);
        p2 = p2vec(t+1);
        u1 = u1vec(t+1);
        u2 = u2vec(t+1);
        v1 = v1vec(t+1);
        v2 = v2vec(t+1);

        p1f = 0;
        p2f = 0;
        v1f = 0;
        v2f = 0;

        for lu = 1:nu

            p1f = p1f + Pu(ju,lu)*intf3(knotsy,knotsy,knotsx,reshape(p1mat(ju,lu,:,:,:),[ny ny nx]),y1,y2,x);
            p2f = p2f + Pu(ju,lu)*intf3(knotsy,knotsy,knotsx,reshape(p2mat(ju,lu,:,:,:),[ny ny nx]),y1,y2,x);
            v1f = v1f + Pu(ju,lu)*intf3(knotsy,knotsy,knotsx,reshape(v1mat(ju,lu,:,:,:),[ny ny nx]),y1,y2,x);
            v2f = v2f + Pu(ju,lu)*intf3(knotsy,knotsy,knotsx,reshape(v2mat(ju,lu,:,:,:),[ny ny nx]),y1,y2,x);

        end

        e1vec(t+1) = -p1 + bet*p1f + kap1*u1 + kap1*((rho+eta)*y1 + (1-gam)*(1-rho)*(y1-y2));
        e2vec(t+1) = -p2 + bet*p2f + kap2*u2 + kap2*((rho+eta)*y2 - gam*(1-rho)*(y1-y2));    
        e3vec(t+1) = -1 - (gam*(eta+rho)*(y1 - (1-gam)/gam/(eta+rho)*u1)^2 ...
            + (1-gam)*(eta+rho)*(y2 + 1/(eta+rho)*u2)^2 ...
            + gam*(1-gam)*(1-rho)*(y1-y2)^2 ...
            + gam*sig/kap1*p1^2 + (1-gam)*sig/kap2*p2^2)/v1 + bet*v1f/v1;
        e4vec(t+1) = -1 - ((1-gam)*(eta+rho)*(y2 - gam/(1-gam)/(eta+rho)*u2)^2 ...
            + gam*(eta+rho)*(y1 + 1/(eta+rho)*u1)^2 ...
            + gam*(1-gam)*(1-rho)*(y2-y1)^2 ...
            + gam*sig/kap1*p1^2 + (1-gam)*sig/kap2*p2^2)/v2 + bet*v2f/v2;

    end

    er_mean(1,i) = mean(abs(e1vec(drop+2:T+1)));
    er_mean(2,i) = mean(abs(e2vec(drop+2:T+1)));
    er_mean(3,i) = mean(abs(e3vec(drop+2:T+1)));
    er_mean(4,i) = mean(abs(e4vec(drop+2:T+1)));
    er_max(1,i) = max(abs(e1vec(drop+2:T+1)));
    er_max(2,i) = max(abs(e2vec(drop+2:T+1)));
    er_max(3,i) = max(abs(e3vec(drop+2:T+1)));
    er_max(4,i) = max(abs(e4vec(drop+2:T+1)));

    sd(1,i) = std(y1vec(drop+2:T+1)-(1-gam)/gam/(eta+rho)*u1vec(drop+2:T+1));
    sd(2,i) = std(y2vec(drop+2:T+1)+1/(eta+rho)*u2vec(drop+2:T+1));
    sd(3,i) = std(y1vec(drop+2:T+1)-y2vec(drop+2:T+1));
    sd(4,i) = std(p1vec(drop+2:T+1));
    sd(5,i) = std(p2vec(drop+2:T+1));
    sd(6,i) = std(y2vec(drop+2:T+1)-gam/(1-gam)/(eta+rho)*u2vec(drop+2:T+1));
    sd(7,i) = std(y1vec(drop+2:T+1)+1/(eta+rho)*u1vec(drop+2:T+1));
    
    ev(1,i) = -(gam*(eta+rho)*sd(1,i)^2 ...
            + (1-gam)*(eta+rho)*sd(2,i)^2 ...
            + gam*(1-gam)*(1-rho)*sd(3,i)^2 ...
            + gam*sig/kap1*sd(4,i)^2 + (1-gam)*sig/kap2*sd(5,i)^2)/(1-bet);
    ev(2,i) = -((1-gam)*(eta+rho)*sd(6,i)^2 ...
            + gam*(eta+rho)*sd(7,i)^2 ...
            + gam*(1-gam)*(1-rho)*sd(3,i)^2 ...
            + gam*sig/kap1*sd(4,i)^2 + (1-gam)*sig/kap2*sd(5,i)^2)/(1-bet);
    
    ivec1(:,i) = ivec(drop+2:T+1);
    xvec1(:,i) = xvec(drop+2:T+1);
    zvec1(:,i) = zvec(drop+2:T+1);
    
end

er1(1,1) = log10(mean(er_mean(1,:)));
er1(2,1) = log10(mean(er_mean(2,:)));
er1(3,1) = log10(mean(er_mean(3,:)));
er1(4,1) = log10(mean(er_mean(4,:)));
er1(1,2) = log10(max(er_max(1,:)));
er1(2,2) = log10(max(er_max(2,:)));
er1(3,2) = log10(max(er_max(3,:)));
er1(4,2) = log10(max(er_max(4,:)));

sd1(1:7) = mean(sd(1:7,:),2);
ev1(1:2) = mean(ev(1:2,:),2);
ev1(3) = v1mat(5,5,midy,midy,midx);
ev1(4) = v2mat(5,5,midy,midy,midx);

disp(er1');
disp(sd1);
disp(ev1);