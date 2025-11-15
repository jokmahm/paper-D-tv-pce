function y = shifted_legendre(n,a,b,x)

if iscolumn(x)
    z = x;
else
    z = x';
end
pl = plgndr(2*((z-a)/(b-a))-1,n);
y = pl(:,end);