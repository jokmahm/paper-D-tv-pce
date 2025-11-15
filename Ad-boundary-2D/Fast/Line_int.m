function [ux_x0, ux_x1, uy_0y, uy_1y, dux0, dux1, du0y, du1y] = Line_int(C, p, Bx_x0, Bx_x1,Bx_0y, Bx_1y, By_x0, By_x1, By_0y, By_1y)
m = length(Bx_x0(1,:));
Np = nchoosek(p+2,2);

ux_x0 = zeros(1,m);
ux_x1 = zeros(1,m);
for i = 1:Np
   ux_x0 = ux_x0+C(i)*Bx_x0(i,:);
   ux_x1 = ux_x1+C(i)*Bx_x1(i,:);
end

uy_0y = zeros(1,m);
uy_1y = zeros(1,m);
for i = 1:Np
   uy_0y = uy_0y+C(i)*By_0y(i,:);
   uy_1y = uy_1y+C(i)*By_1y(i,:);
end

uy_x0 = zeros(1,m);
uy_x1 = zeros(1,m);
for i = 1:Np
   uy_x0 = uy_x0+C(i)*By_x0(i,:);
   uy_x1 = uy_x1+C(i)*By_x1(i,:);
end

ux_0y = zeros(1,m);
ux_1y = zeros(1,m);
for i = 1:Np
   ux_0y = ux_0y+C(i)*Bx_0y(i,:);
   ux_1y = ux_1y+C(i)*Bx_1y(i,:);
end

h = 1e-12;
dux0 = 1./sqrt(ux_x0.^2+uy_x0.^2+h);
dux1 = 1./sqrt(ux_x1.^2+uy_x1.^2+h);
du0y = 1./sqrt(ux_0y.^2+uy_0y.^2+h);
du1y = 1./sqrt(ux_1y.^2+uy_1y.^2+h);