function I = Gauss_Quadrature1D(f,a,b,k)
% 2 D integration with Gauss Quaderature rule
% for x in [a c] and y in [c d]

% points and weigths
[p,w] = grule(k);    

% Coordinate transformation
xp(1:k) = 0.5*((b-a).*p(1:k)+(b+a));

% Evaluation
I = (b-a)/2*sum(f(xp).*w);