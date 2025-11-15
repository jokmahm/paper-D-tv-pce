function y = Test_1D(z)
n = length(z);
y = zeros(1,n);
for i = 1:n
    if z(i)<=0.5 
        y(i) = 1; % z(i)-1; 
    else
        y(i) = 0; % z(i);
    end
end