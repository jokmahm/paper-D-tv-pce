function gamma = Mass(p,l,u)
global idx
Np = nchoosek(p+2,2);
gamma = zeros(Np,1);
for i = 1:Np
    gamma(i) =  ((u-l)./(2*idx(i,1)+1))*((u-l)./(2*idx(i,2)+1));
end