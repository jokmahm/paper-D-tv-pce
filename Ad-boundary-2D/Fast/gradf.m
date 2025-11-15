function [fx_0y, fx_1y, fy_x0, fy_x1, gx0, gx1, g0y, g1y] = gradf(f1,f2,N,a,b,c,d)
[xp,~] = grule(N);
xP(1:N) = 0.5*((b-a).*xp(1:N)+(b+a));
yP(1:N) = 0.5*((d-c).*xp(1:N)+(d+c));
h = 1e-5;

fx_x0 = (f1(xP+h,c)-f1(xP,c))/h;
fx_x1 = (f1(xP+h,d)-f1(xP,d))/h;

fy_0y = (f2(a,yP+h)-f2(a,yP))/h;
fy_1y = (f2(b,yP+h)-f2(b,yP))/h;

fy_x0 = (f1(xP,c+h)-f1(xP,c))/h;
fy_x1 = (f1(xP,d+h)-f1(xP,d))/h;

fx_0y = (f2(a+h,yP)-f2(a,yP))/h;
fx_1y = (f2(b+h,yP)-f2(b,yP))/h;

gx0 = sqrt(fx_x0.^2+fy_x0.^2);
gx1 = sqrt(fx_x1.^2+fy_x1.^2);
g0y = sqrt(fx_0y.^2+fy_0y.^2);
g1y = sqrt(fx_1y.^2+fy_1y.^2);