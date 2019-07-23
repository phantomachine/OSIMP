% clear all;

ny = length(knotsy);

T = 21;
y1vec = zeros(T+1,1);
y2vec = zeros(T+1,1);
p1vec = zeros(T+1,1);
p2vec = zeros(T+1,1);
v1vec = zeros(T+1,1);
v2vec = zeros(T+1,1);
ivec = 5*ones(T+1,1);
ivec(2) = 8; % positive markup shock

for i = 1:T
    
    ju = ivec(i+1);
    
    y1vec(i+1) = intf2(knotsy,knotsy,reshape(y1mat(ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    y2vec(i+1) = intf2(knotsy,knotsy,reshape(y2mat(ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    p1vec(i+1) = intf2(knotsy,knotsy,reshape(p1mat(ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    p2vec(i+1) = intf2(knotsy,knotsy,reshape(p2mat(ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    v1vec(i+1) = intf2(knotsy,knotsy,reshape(v1mat(ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    v2vec(i+1) = intf2(knotsy,knotsy,reshape(v2mat(ju,:,:),[ny ny]),y1vec(i),y2vec(i));

end

y1vec = y1vec(2:end);
y2vec = y2vec(2:end);
p1vec = p1vec(2:end);
p2vec = p2vec(2:end);
v1vec = v1vec(2:end);
v2vec = v2vec(2:end);
ttvec = y1vec-y2vec;

% figure;
% subplot(321);
% plot([1:T],y1vec(1:end));
% xlim([1 T]);
% subplot(322);
% plot([1:T],y2vec(1:end));
% xlim([1 T]);
% subplot(323);
% plot([1:T],p1vec(1:end));
% xlim([1 T]);
% %ylim([-0.01 0.02]);
% subplot(324);
% plot([1:T],p2vec(1:end));
% xlim([1 T]);
% subplot(325);
% plot([1:T],ttvec(1:end));
% xlim([1 T]);
% 
% eval(['save ./ifsc_k3_ny5nx5rho1.0_neu', num2str(neu,'%1.2f'), '.mat y1vec y2vec p1vec p2vec v1vec v2vec ttvec;']);