function gamma = Mass(p,l1,u1,l2,u2)
global idx
Np = nchoosek(p+2,2);
gamma = zeros(Np,1);
for i = 1:Np
    gamma(i) =  ((u1-l1)./(2*idx(i,1)+1))*((u2-l2)./(2*idx(i,2)+1));
end