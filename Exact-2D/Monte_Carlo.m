function [mu,var] = Monte_Carlo(f,N)

y= zeros(N,1);
for i = 1:N
    z1 = rand;
    z2 = rand;
    y(i) = f(z1,z2);    
end
mu = mean(y);

v = zeros(N,1);
for i = 1:N
    v(i) = (y(i)-mu)^2;
end
var = sum(v)/(N-1);