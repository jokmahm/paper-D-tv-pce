function I = Gauss_Quadrature2D(f,a,b,c,d,k)
% 2 D integration with Gauss Quaderature rule
% for x in [a c] and y in [c d]

[p,w] = grule(k);    % points and weigths

% Coordinate transformation
xp(1:k) = 0.5*((b-a).*p(1:k)+(b+a));
yp(1:k) = 0.5*((d-c).*p(1:k)+(d+c));

% Evaluation

% Vectorized
[Xp,Yp] = meshgrid(xp,yp);
W = w'*w;
E = f(Xp,Yp).*W;
I = (b-a)*(d-c)/4*sum(sum(E));

% Pointwise
% E = zeros(k,k);
% for i = 1:k
%     for j = 1:k
%         E(i,j) = f(xp(i),yp(j))*w(i)*w(j);
%     end
% end
% 
% I = (b-a)*(d-c)/4*sum(sum(E));