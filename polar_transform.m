function [u,G_base] = polar_transform(N,X)
% polarization transformation 
% input: N 
% ouput: u transofrmed 
n = log2(N); 
G_2 = [1 0; 1 1]; 
G_base = G_2; 

for i=1:(n-1)
G_base = kron(G_base,G_2); 
end 

u = mod(X*G_base,2); 

end