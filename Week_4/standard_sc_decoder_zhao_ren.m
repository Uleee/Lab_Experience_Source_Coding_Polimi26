function [u,u_hat,idx]=standard_sc_decoder_zhao_ren(P,N,u_hat,idx,h_thr,u_G)

if N==1 
 % stage 1 entropy calculation, ML descicsion
p_0  = P(1,1);  
p_1 = P(1,2); 
if p_0<=0 || p_0>=1 
     h=0;
else
    h = -p_0*log2(p_0) - (p_1)*log2(p_1); 
end


if h>=h_thr
    u = u_G(idx); 
    idx = idx + 1;
else
    if p_0 >= p_1
    u = 0;
    else 
    u = 1; 
    end 
end

u_hat = [u_hat, u];

else 

N_top = N/2; 
P_top = P(1:N/2,:); 
P_bottom = P( N/2+1:end,:); 

% minus ch
P_left = zeros(N_top, 2);
for k = 1:N_top 
 P_left(k, 1) = 1 * (P_top(k, 1)*P_bottom(k, 1) + P_top(k, 2)*P_bottom(k, 2));
 P_left(k, 2) = 1 * (P_top(k, 2)*P_bottom(k, 1) + P_top(k, 1)*P_bottom(k, 2));
 sum_P_left = P_left(k, 1) + P_left(k, 2);
        if sum_P_left > 0
            P_left(k, :) = P_left(k, :) / sum_P_left;
        end
end

[upper_u, u_hat,idx] = standard_sc_decoder_zhao_ren(P_left,N_top,u_hat,idx,h_thr,u_G); 

% plus channel
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
[lower_u, u_hat,idx] = standard_sc_decoder_zhao_ren(P_right,N_top,u_hat,idx,h_thr,u_G); 

u=[xor(upper_u,lower_u), lower_u]; 
end 

