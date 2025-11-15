clc; clear;
% 1D minimization of Total Variation (1D-TV)
% Exact form of Legender construction

l = 0;                           % parameter uncertainty lower bound
u = 1;                           % parameter uncertainty upper bound
p = 5;                          % degree of the polynomial
N = p+10;                        % order Gauss quaderature integration rule
f = @(x) x.^5-4*x.^3+2; %sin(x);
df = @(x) 5*x.^4-12*x.^2; %cos(x);
ddf = @(x) 20*x.^3-24*x; %-sin(x);
b = zeros(p+1,1);                % right hand side
[Bp, ~, ~, W, V] = Guass_quad_integration(p,l,u,N,f);
for i = 1:p+1
   I = Bp(:,i).*V'.*W';
   b(i) = (u-l)/2*sum(I);
end 
gamma = (u-l)./(2*(0:p)+1);       % orthogonality coefficients
M = diag(gamma);                  % Mass matrix
Xi = M\b;  

y = linspace(l,u,100);

F = f(y');
dF = df(y');
ddF = ddf(y');

B = basis(p,l,u,y);
B_x = basis_x(p,l,u,y);
B_xx = basis_xx(p,l,u,y);
g = B*Xi;
dg = B_x*Xi; 
ddg = B_xx*Xi;

figure
plot(y,F)
hold on
plot(y,g,'r--')
legend('exact','approx')
title('function f')
fprintf('Function Error: %0.4e \n', norm(F-g))
figure
plot(y,dF)
hold on
plot(y,dg,'r--')
legend('exact','approx')
title('function df')
fprintf('First drivative Error: %0.4e \n', norm(dF-dg))
figure
plot(y,ddF)
hold on
plot(y,ddg,'r--')
legend('exact','approx')
title('function ddf')
fprintf('Second drivative Error: %0.4e \n', norm(ddF-ddg))