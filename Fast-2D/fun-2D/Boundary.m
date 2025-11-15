function y = Boundary(a)
f = @ (x) (1-x).^(2/3).*exp(-x);
y = 1/integral(f,0,1);