clear; clc;
global idx
l = 0;                           % parameter uncertainty lower bound
u = 1;                           % parameter uncertainty upper bound
p = 16;                           % degree of the polynomial
n = 20;                          % number of discritization
x1 = linspace(l,u,n);            % discritization for simpson integration
x2 = linspace(l,u,n);
[X1, X2] = meshgrid(x1,x2);
% f = @(x,y) sin(x+y);
f = @(x,y) Mc_2D(0.5,0.5,x,y,3);
Z = f(X1,X2);
Np = nchoosek(p+2,2);
b = zeros(Np,1);                  % right hand side

idx = zeros(Np,2);
k = 1;
for i = 0:p
    for j = 0:p
        if(i+j<=p)
            idx(k,:) = [i,j];
            k = k+1;
        end
    end
end

for i = 1:Np
   g = @(x,y) Phi(i,l,u,x,y);
   h = @(x,y) g(x,y) .* f(x,y);
   b(i) = integral2(h,l,u,l,u);
   % b(i) = Gauss_Quadrature2D(h,l,u,l,u,10);
end 

gamma = Mass(p,l,u);              % orthogonality coefficients
M = diag(gamma);                  % Mass matrix
c = M\b;                          % coefficients

Q = zeros(size(Z));
B = basis(p,l,u,X1,X2);
for i = 1:Np
    Q = Q+c(i)*B(:,:,i);
end

subplot(1,3,1), surf(X1,X2,Z)
xlabel('x')
ylabel('y')
zlabel('Exact')

subplot(1,3,2), surf(X1,X2,Q)
xlabel('x')
ylabel('y')
zlabel('Approximate')

subplot(1,3,3), surf(X1,X2,abs(Z-Q))
xlabel('x')
ylabel('y')
zlabel('Error')