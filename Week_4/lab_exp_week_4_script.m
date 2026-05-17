% lab experience: [2025] Entropy Polarization-Based Data Compression Without Frozen Set Construction
clc
clear 
% generate brv 
% apply polarization transformation 
% 
%% Task 1 
clc
clear 
%
N = 1024; 
p = 0.2; % Bernoulli prob

X = mod(rand(1,N) < p,2); 
% probability of X != true bit is 0.2, X == true bit = 1-0.2 = 0.8;
[u_true,G_base] = polar_transform(N,X); 
y = zeros(1,N); 
% LLR = zeros(1,N)+log((1-p)/p);
% frozen_ind = zeros(1, N);
P = zeros(N,2); 
P(:, 1) = 1-p; 
P(:, 2) = p; 

idx=0; 
[~,~,h_set,~,T]=standard_sc_encoder_zhao_ren(P,N,[],[],u_true,idx,[]);

if isempty(T) % sometimes T can be empty (for allzero u_true)
    h_thr = 0; 
else
    h_thr = min(h_set(T)); 
end

G = find(h_set >= h_thr);
u_G = u_true(G); 

idx = 1; % tau pointter

[u,u_hat,idx]=standard_sc_decoder_zhao_ren(P,N,[],idx,h_thr,u_G); 

X_est = mod(u_hat*G_base,2);

check=(X_est==X);

if all(check)
    disp('The sequences match perfectly');
else
    disp('There are some errors in the reconstruction.');
end

%% E{R} calculation 


clc
clear
prob = 0.001:0.01:0.5; 
n=10; 
N = 2^n; 
ER_ren_zhao = zeros(1, numel(prob));

G_2 = [1 0; 1 1]; 
G_base = G_2; 
for i=1:(n-1)
G_base = kron(G_base,G_2); 
end 

for k = 1:numel(prob)
p = prob(k);


for m = 1:1000
    
X = (rand(1,N) < p);

u_true = mod(X*G_base,2);    

idx = 0; 

P = zeros(N,2); 
P(:, 1) = 1-p; 
P(:, 2) = p; 

[~,~,h_set,~,T]=standard_sc_encoder_zhao_ren(P,N,[],[],u_true,idx,[]);

if isempty(T) % sometimes T can be empty (for allzero u_true)
    h_thr = 0; 
else
    h_thr = min(h_set(T)); 
end

G = find(h_set >= h_thr);
u_G = u_true(G); 

idx_tau = 1; % tau pointter

[u,u_hat,idx]=standard_sc_decoder_zhao_ren(P,N,[],idx_tau,h_thr,u_G); 

EL(m) = numel(u_G); 

end

EL = mean(EL); 
ER_ren_zhao(1,k) = EL/N; 
end 


h2p = -prob.*log2(prob)-(1-prob).*log2(1-prob); 

%% 
plot(prob,ER_ren_zhao,'r')
hold on 
plot(prob,h2p,'b')
hold on 
plot(prob,ER_korada,'g')




%% NON STANDARD SC ENCODER/DECODER 

clc
clear 
%
N = 1024; 
p = 0.2; % Bernoulli prob

X = mod(rand(1,N) < p,2); 
% probability of X != true bit is 0.2, X == true bit = 1-0.2 = 0.8;
[u_true,G_base] = polar_transform(N,X); 
y = zeros(1,N); 
LLR = zeros(1,N)+log((1-p)/p);

idx = 0; 

[~,~,T,~,LLRset] = non_standard_sc_encoder_zhao_ren(LLR,N,[],u_true,idx,[],[]);

if isempty(T)  
    llr_thr = 0; 
else
    llr_thr = max(abs(LLRset(T))); 
end

G = find(abs(LLRset) <= llr_thr);
u_G = u_true(G); 

idx_tau = 1; 
[u,u_hat,idx,LLRset2] = non_standard_sc_decoder_zhao_ren(LLR,N,[],idx_tau,llr_thr,u_G,[]); 

X_est = mod(u_hat*G_base,2);

check=(X_est==X);

if all(check)
    disp('The sequences match perfectly');
else
    disp('There are some errors in the reconstruction.');
end


%% E[R] calculation 
clc
clear
prob = 0.001:0.01:0.5; 
n=10; 
N = 2^n; 
ER_ren_zhao_soft = zeros(1, numel(prob));

G_2 = [1 0; 1 1]; 
G_base = G_2; 
for i=1:(n-1)
G_base = kron(G_base,G_2); 
end 

for k = 1:numel(prob)
p = prob(k);


for m = 1:1000
    
X = (rand(1,N) < p);

u_true = mod(X*G_base,2);    
LLR = zeros(1,N)+log((1-p)/p);
idx = 0; 

[~,~,T,~,LLRset] = non_standard_sc_encoder_zhao_ren(LLR,N,[],u_true,idx,[],[]);

if isempty(T)  
    llr_thr = 0; 
else
    llr_thr = max(abs(LLRset(T))); 
end

G = find(abs(LLRset) <= llr_thr);
u_G = u_true(G); 

idx_tau = 1; 
[u,u_hat,idx,LLRset2] = non_standard_sc_decoder_zhao_ren(LLR,N,[],idx_tau,llr_thr,u_G,[]); 

EL(m) = numel(u_G); 

end

EL = mean(EL); 
ER_ren_zhao_soft(1,k) = EL/N; 
end 

h2p = -prob.*log2(prob)-(1-prob).*log2(1-prob); 

%% 
plot(prob,ER_ren_zhao_soft,'r')
hold on 
plot(prob,h2p,'b')





%% Section III. D fixed threshold Let's implement it with LLRs 
clc
clear 
N = 256; 
p = 0.2; 
X = mod(rand(1,N)<p,2); 
[u_true,G_base] = polar_transform(N,X); 
y=zeros(1,N); 
LLR = zeros(1,N)+log((1-p)/p);

idx = 0; 

[~,~,T,~,LLRset] = non_standard_sc_encoder_zhao_ren(LLR,N,[],u_true,idx,[],[]);

eps_fixed = 1/(log2(N)); 
fixed_llr_thr = log((1-eps_fixed)/(eps_fixed)); 

if isempty(T)  
    llr_thr = 0; 
else
    llr_thr = max(abs(LLRset(T))); 
end

G = find(abs(LLRset) <= llr_thr);
G_fix = find(abs(LLRset) <= fixed_llr_thr);
G_fixed_star = intersect(setdiff(G, G_fix), T);
u_G_fix=u_true(G_fix); 

idx = 0; 
idx_u = 0;
[u,u_hat,idx] = non_standard_sc_decoder_zhao_ren_fixed_thr(LLR,N,[],idx,idx_u,llr_thr,fixed_llr_thr,u_G_fix,G_fixed_star);

X_est = mod(u_hat*G_base,2);

check=(X_est==X);

if all(check)
    disp('The sequences match perfectly');
else
    disp('There are some errors in the reconstruction.');
end


%% E[R] calculation 
clc
clear
prob = 0.001:0.01:0.5; 
n=10; 
N = 2^n; 
ER_ren_zhao_soft_fixed = zeros(1, numel(prob));

G_2 = [1 0; 1 1]; 
G_base = G_2; 
for i=1:(n-1)
G_base = kron(G_base,G_2); 
end 

for k = 1:numel(prob)
p = prob(k);


for m = 1:1000
    
X = (rand(1,N) < p);

u_true = mod(X*G_base,2);    
LLR = zeros(1,N)+log((1-p)/p);
idx = 0; 
[~,~,T,~,LLRset] = non_standard_sc_encoder_zhao_ren(LLR,N,[],u_true,idx,[],[]);

eps_fixed = 1/(log2(N)); 
fixed_llr_thr = log((1-eps_fixed)/(eps_fixed)); 

if isempty(T)  
    llr_thr = 0; 
else
    llr_thr = max(abs(LLRset(T))); 
end

G = find(abs(LLRset) <= llr_thr);
G_fix = find(abs(LLRset) <= fixed_llr_thr);
G_fixed_star = intersect(setdiff(G, G_fix), T);
u_G_fix=u_true(G_fix); 

idx = 0; 
idx_u = 0;
[u,u_hat,idx] = non_standard_sc_decoder_zhao_ren_fixed_thr(LLR,N,[],idx,idx_u,llr_thr,fixed_llr_thr,u_G_fix,G_fixed_star);

EL(m) = log2(N)+numel(G_fix)+numel(G_fixed_star)*(log2(N));  

end

EL = mean(EL); 
ER_ren_zhao_soft_fixed(1,k) = EL/N; 
end 

h2p = -prob.*log2(prob)-(1-prob).*log2(1-prob); 

%% 
plot(prob,ER_ren_zhao,'r')
hold on 
plot(prob,ER_ren_zhao_soft,'b')
hold on 
plot(prob,ER_ren_zhao_soft_fixed,'g')
hold on 
plot(prob,h2p,'p')

















