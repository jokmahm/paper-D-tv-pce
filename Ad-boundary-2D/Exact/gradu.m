function [u_x, u_y] = gradu(C,p,B_x,B_y)

[n,m] = size(B_x(:,:,1));
u_x = zeros(n,m);
u_y = zeros(n,m);

Np = nchoosek(p+2,2);
for i = 1:Np
   u_x = u_x+C(i)*B_x(:,:,i);
   u_y = u_y+C(i)*B_y(:,:,i);
end