knotsy = knots1;
ny = length(knotsy);
nx = length(knotsx);

T = 21;
y1vec = zeros(T+1,1); % y1(t-1)
y2vec = zeros(T+1,1); % y2(t-1)
xvec = .5*ones(T+1,1); % y2(t-1)
zvec = ones(T+1,1); % y2(t-1)
p1vec = zeros(T+1,1); % p1(t)
p2vec = zeros(T+1,1);
v1vec = zeros(T+1,1);
v2vec = zeros(T+1,1);
ivec = 5*ones(T+1,1);
ivec(2) = 8; % positive markup shock

for i = 1:T
    
    ku = ivec(i);
    ju = ivec(i+1);
    
    y1vec(i+1) = intf3(knotsy,knotsy,knotsx,reshape(y1mat(ku,ju,:,:,:),[ny ny nx]),y1vec(i),y2vec(i),xvec(i));
    y2vec(i+1) = intf3(knotsy,knotsy,knotsx,reshape(y2mat(ku,ju,:,:,:),[ny ny nx]),y1vec(i),y2vec(i),xvec(i));
    p1vec(i+1) = intf3(knotsy,knotsy,knotsx,reshape(p1mat(ku,ju,:,:,:),[ny ny nx]),y1vec(i),y2vec(i),xvec(i));
    p2vec(i+1) = intf3(knotsy,knotsy,knotsx,reshape(p2mat(ku,ju,:,:,:),[ny ny nx]),y1vec(i),y2vec(i),xvec(i));    
    v1vec(i+1) = intf3(knotsy,knotsy,knotsx,reshape(v1mat(ku,ju,:,:,:),[ny ny nx]),y1vec(i),y2vec(i),xvec(i));
    v2vec(i+1) = intf3(knotsy,knotsy,knotsx,reshape(v2mat(ku,ju,:,:,:),[ny ny nx]),y1vec(i),y2vec(i),xvec(i));    
    zvec(i+1)  = intf3(knotsy,knotsy,knotsx,reshape(zmat(ku,ju,:,:,:),[ny ny nx]),y1vec(i),y2vec(i),xvec(i));    
    xvec(i+1)  = intf3(knotsy,knotsy,knotsx,reshape(xmat(ku,ju,:,:,:),[ny ny nx]),y1vec(i),y2vec(i),xvec(i));    
    
    u1vec(i+1) = Gu(ivec(i),1);
    u2vec(i+1) = Gu(ivec(i),2);

end

y1vec = y1vec(2:end);
y2vec = y2vec(2:end);
p1vec = p1vec(2:end);
p2vec = p2vec(2:end);
v1vec = v1vec(2:end);
v2vec = v2vec(2:end);
ttvec = y1vec-y2vec;

figure;
subplot(321);
plot([1:T],y1vec(1:end),'b-');
hold on;
plot([1:T],y2vec(1:end),'r--');
xlim([1 T]);
subplot(322);
plot([1:T],v1vec(1:end));
xlim([1 T]);
subplot(323);
plot([1:T],p1vec(1:end),'b-');
hold on;
plot([1:T],p2vec(1:end),'r--');
xlim([1 T]);
subplot(324);
plot([1:T],v2vec(1:end));
xlim([1 T]);
subplot(325);
plot([1:T],ttvec(1:end));
xlim([1 T]);