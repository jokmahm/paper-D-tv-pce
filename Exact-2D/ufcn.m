function z = ufcn(p,C,B_x,B_y)
[u_x, u_y] = gradu(C,p,B_x,B_y);
z = 1./sqrt(u_x.^2+u_y.^2+1e-12);