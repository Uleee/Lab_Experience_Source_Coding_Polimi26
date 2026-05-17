function [u,u_hat,idx,LLRset] = non_standard_sc_decoder_zhao_ren(LLR,N,u_hat,idx,llr_thr,u_G,LLRset)
if N==1 

    if abs(LLR)<=llr_thr
    u = u_G(idx); 
    idx = idx + 1;
else
    if LLR >= 0
    u = 0;
    else 
    u = 1; 
    end 
    end
    
u_hat = [u_hat,u];
LLRset = [LLRset, LLR];

else 
    
N_top = N/2; 

LLR_top = LLR(1:N/2); 
LLR_bottom = LLR(N/2+1:end); 

LLR_d_minus = 2*atanh( tanh(LLR_top/2) .* tanh(LLR_bottom/2) );

[upper_u, u_hat,idx,LLRset] = non_standard_sc_decoder_zhao_ren(LLR_d_minus,N_top,u_hat,idx,llr_thr,u_G,LLRset); 

LLR_d_plus =  (-1).^(upper_u) .* LLR_top + LLR_bottom;

[lower_u,u_hat, idx,LLRset] = non_standard_sc_decoder_zhao_ren(LLR_d_plus,N_top,u_hat,idx,llr_thr,u_G,LLRset); 

u=[xor(upper_u,lower_u), lower_u]; 

end 
