function [N_St] = Normalizefunction(M)

%Normalize Between '0 to 1'
[m,n]=size(M);
x=int64(M);
M=double(x);

for i=1:m
    for j=1:n
        N_St(i,j)=0.00392*M(i,j);
        
    end
end
end

