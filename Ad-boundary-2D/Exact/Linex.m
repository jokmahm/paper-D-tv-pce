function [Bx0, Bx1, B0y, B1y, w] = Linex(p,a,b,c,d,k)
[P,w] = grule(k);    % points and weigths
Np = nchoosek(p+2,2);

xp(1:k) = 0.5*((b-a).*P(1:k)+(b+a));
yp(1:k) = 0.5*((d-c).*P(1:k)+(d+c));

Bx0 = zeros(Np,k);
Bx1 = zeros(Np,k);
B0y = zeros(Np,k);
B1y = zeros(Np,k);

for i = 1:Np
    Bx0(i,:) = Phi(i-1,a,b,c,d,xp,c);
    Bx1(i,:) = Phi(i-1,a,b,c,d,xp,d);
    
    B0y(i,:) = Phi(i-1,a,b,c,d,a,yp);
    B1y(i,:) = Phi(i-1,a,b,c,d,b,yp);
end