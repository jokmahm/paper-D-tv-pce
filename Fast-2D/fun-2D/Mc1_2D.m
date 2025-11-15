function y = Mc1_2D(a,t,z1,z2)
[m,n] = size(z1);
y = zeros(m,n);
for i = 1:m
    for j = 1:n
        if z2(i,j)<=0.5 
            y(i,j) = (z1(i,j)/2)*exp(t)*(1-a)^(2/3)*exp(-a);
        else
            y(i,j) = z1(i,j)*exp(t)*(1-a)^(2/3)*exp(-a);
        end
    end
end