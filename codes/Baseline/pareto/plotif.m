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
    
    ku = ivec(i);
    ju = ivec(i+1);
    
    y1vec(i+1) = intf2(knotsy,knotsy,reshape(y1mat(ku,ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    y2vec(i+1) = intf2(knotsy,knotsy,reshape(y2mat(ku,ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    p1vec(i+1) = intf2(knotsy,knotsy,reshape(p1mat(ku,ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    p2vec(i+1) = intf2(knotsy,knotsy,reshape(p2mat(ku,ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    v1vec(i+1) = intf2(knotsy,knotsy,reshape(v1mat(ku,ju,:,:),[ny ny]),y1vec(i),y2vec(i));
    v2vec(i+1) = intf2(knotsy,knotsy,reshape(v2mat(ku,ju,:,:),[ny ny]),y1vec(i),y2vec(i));

end

y1vec = y1vec(2:end);
y2vec = y2vec(2:end);
p1vec = p1vec(2:end);
p2vec = p2vec(2:end);
v1vec = v1vec(2:end);
v2vec = v2vec(2:end);
ttvec = y1vec-y2vec;