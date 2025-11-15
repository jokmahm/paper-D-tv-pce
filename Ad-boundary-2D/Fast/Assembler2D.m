function J = Assembler2D(p, Xi, B_x, B_y, W, a, b, c, d)
du = ufcn(p,Xi,B_x,B_y);                 % a(Xi)
Np = nchoosek(p+2,2);
J = zeros(Np,Np);
for i = 1:Np
    for j = 1:Np
        dy1 = B_x(:,:,i).*B_x(:,:,j);
        dy2 = B_y(:,:,i).*B_y(:,:,j);
        h = du.*(dy1+dy2);
        E = h.*W;
        J(i,j) = ((b-a)/2)*((d-c)/2)*sum(sum(E));   % gauss integration rule
    end
end