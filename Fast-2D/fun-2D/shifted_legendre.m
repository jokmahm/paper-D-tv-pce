function y = shifted_legendre(n,a,b,x)

%[y, ~] = legendrePoly(n,2*((x-a)/(b-a))-1);
[r,s] = size(x);
y = zeros(r,s);
for i = 1:s
    h = plgndr(2*((x(:,i)-a)/(b-a))-1,n);
    y(:,i) = h(:,end);
end