function [u,u_hat,T,idx,LLRset] = non_standard_sc_encoder_zhao_ren(LLR,N,u_hat,u_true,idx,T,LLRset)

if N==1
    idx = idx+1; 

if LLR >= 0
   u_ml = 0;
else
   u_ml = 1;
end
        
if u_true(idx) ~= u_ml
T = [T, idx]; 
end
        
u = u_true(idx); 
u_hat = [u_hat, u]; 
LLRset = [LLRset, LLR];

else 
    
N_top = N/2; 

LLR_top = LLR(1:N/2); 
LLR_bottom = LLR(N/2+1:end); 

LLR_d_minus = 2*atanh( tanh(LLR_top/2) .* tanh(LLR_bottom/2) );

[upper_u, u_hat,T,idx,LLRset] = non_standard_sc_encoder_zhao_ren(LLR_d_minus,N_top,u_hat,u_true,idx,T,LLRset); 

LLR_d_plus =  (-1).^(upper_u) .* LLR_top + LLR_bottom;

[lower_u,u_hat, T,idx,LLRset] = non_standard_sc_encoder_zhao_ren(LLR_d_plus,N_top,u_hat,u_true,idx,T,LLRset); 

u=[xor(upper_u,lower_u), lower_u]; 

end 

