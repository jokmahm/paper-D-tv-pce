function y = Test_2D(z1,z2)
[m,n] = size(z1);
y = zeros(m,n);
for i = 1:m
    for j = 1:n
        if z2(i,j)<=0.5 % && z2(i,j)>=0.4
              y(i,j) = 0; % z1(i,j)
         else
               y(i,j) = 1; % 2*z1(i,j)
         end
     end
end