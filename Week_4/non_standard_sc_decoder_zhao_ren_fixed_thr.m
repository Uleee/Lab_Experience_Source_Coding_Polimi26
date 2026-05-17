function [u,u_hat,idx,idx_u] = non_standard_sc_decoder_zhao_ren_fixed_thr(LLR,N,u_hat,idx,idx_u,llr_thr,fixed_llr,u_G,G_fix_star)
if N == 1 
idx = idx + 1;

if abs(LLR) <= fixed_llr
idx_u = idx_u + 1;
u = u_G(idx_u); 
else
if LLR >= 0
u = 0;
else 
u = 1; 
end 

if ismember(idx, G_fix_star)
  u = mod(u+1,2);  
end
end
        
u_hat = [u_hat u];
else 
    
N_top = N/2; 

LLR_top = LLR(1:N/2); 
LLR_bottom = LLR(N/2+1:end); 

LLR_d_minus = 2*atanh( tanh(LLR_top/2) .* tanh(LLR_bottom/2) );

[upper_u, u_hat,idx,idx_u] = non_standard_sc_decoder_zhao_ren_fixed_thr(LLR_d_minus,N_top,u_hat,idx,idx_u,llr_thr,fixed_llr,u_G,G_fix_star); 

LLR_d_plus =  (-1).^(upper_u) .* LLR_top + LLR_bottom;

[lower_u,u_hat, idx,idx_u] = non_standard_sc_decoder_zhao_ren_fixed_thr(LLR_d_plus,N_top,u_hat,idx,idx_u,llr_thr,fixed_llr,u_G,G_fix_star); 

u=[xor(upper_u,lower_u), lower_u]; 

end 
