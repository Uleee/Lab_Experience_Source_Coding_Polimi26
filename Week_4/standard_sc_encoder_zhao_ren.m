function [u,u_hat,h_set,idx,T]=standard_sc_encoder_zhao_ren(P,N,h_set,u_hat,u_true,idx,T)
if N==1 
    % stage 1 entropy calculation
idx = idx + 1; 
U_true = u_true(idx); 

p_0  = P(1,1);  
p_1 =  P(1,2); 


if p_0<=0 || p_0>=1 
     h=0;
else
    h = -p_0*log2(p_0) - (p_1)*log2(p_1); 
end

% stage 2 ML descision 

if p_0 >= p_1
   u_ml = 0;
    else 
    u_ml = 1; 
end 

if u_ml~=U_true 
   T = [T idx]; 
end 

u = U_true; 

h_set = [h_set, h]; 
u_hat = [u_hat, u];

else 

N_top = N/2; 
P_top = P(1:N/2,:); 
P_bottom = P( N/2+1:end,:); 

P_left = zeros(N_top, 2);

for k = 1:N_top
 P_left(k, 1) = 1 * (P_top(k, 1)*P_bottom(k, 1) + P_top(k, 2)*P_bottom(k, 2));
 P_left(k, 2) = 1* (P_top(k, 2)*P_bottom(k, 1) + P_top(k, 1)*P_bottom(k, 2));
 sum_P_left = P_left(k, 1) + P_left(k, 2);
 if sum_P_left > 0
    P_left(k, :) = P_left(k, :) / sum_P_left; % normalize to get meaningful numbers 
 end
end

[upper_u, u_hat, h_set,idx,T] = standard_sc_encoder_zhao_ren(P_left,N_top,h_set,u_hat,u_true,idx,T); 


P_right = zeros(N_top, 2);
for k = 1:N_top
if upper_u(k) == 0
  P_right(k, 1) = P_top(k, 1) * P_bottom(k, 1);
  P_right(k, 2) = P_top(k, 2) * P_bottom(k, 2);
else
  P_right(k, 1) = P_top(k, 2) * P_bottom(k, 1);
  P_right(k, 2) = P_top(k, 1) * P_bottom(k, 2);
end
sum_P = P_right(k, 1) + P_right(k, 2);
if sum_P > 0
  P_right(k, :) = P_right(k, :) / sum_P;
end

end
[lower_u, u_hat, h_set,idx,T] = standard_sc_encoder_zhao_ren(P_right,N_top,h_set,u_hat,u_true,idx,T); 

u=[xor(upper_u,lower_u), lower_u]; 
end 


