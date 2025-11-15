function y = Mc(a,t,z)
nz = length(z);
y = zeros(1,nz);
for j = 1:nz
    if z(j)>=0.5
        y(j) = (1/10)*exp(t)*(1-a)^(2/3)*exp(-a);
        %y(j) = (1+1*exp(z(j)))*exp(t)*(1-a)^(2/3)*exp(-a);
    else
        y(j) = exp(t)*(1-a)^(2/3)*exp(-a);
        %y(j) = (1+tan(z(j)).^2)*exp(t)*(1-a)^(2/3)*exp(-a);
    end
end