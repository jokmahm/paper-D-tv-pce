clc; clear;
% 1D minimization of Total Variation (1D-TV)
% Numeric form of Legender construction

p = 99;                         % degree of the polynomial
l = 0;                           % parameter uncertainty lower bound  
u = 1;                           % parameter uncertainty upper bound
lambda = 10*p;                       
N = p+10;                        % order Gauss quaderature integration rule
        
% f = @(x) sin(x);
% f = @(z) Test_1D(z);
f = @(z) Mc(0.5,0.5,z);

b = zeros(p+1,1);                % right hand side
[Bp, Bp_x, W, V] = Guass_quad_integration(p,l,u,N,f);

for i = 1:p+1
   I = Bp(:,i).*V'.*W';
   b(i) = (u-l)/2*sum(I);
end 

gamma = (u-l)./(2*(0:p)+1);       % orthogonality coefficients
M = diag(gamma);                  % Mass matrix
Xi = M\b;                         % coefficients
C = Xi;
r = lambda*b;   % +basis(p,l,u,u)'-basis(p,l,u,l)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = linspace(l,u,200);
B = basis(p,l,u,y);
F = f(y);
U_0 = B*C;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard mean and variance report
mean_pce = C(1);
H = gamma(2:end);
var_pce = H*C(2:end).^2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iter = 10;
Energy = zeros(1,iter);
for k=1:iter                       % non-linear loop
    fprintf('Iteration %d \n', k)         % f', k, norm(r))
    K = Assembeler(Xi,Bp_x,W,p,l,u);
    du = Bp_x*Xi; 
    du(du==0) = 1e-8;
    Energy(k) = norm(du,1)+0.5*lambda*norm(F-U_0,2);
    J = K+lambda*M;
    Xi = J\r;                      % solve for correction
end

figure
plot(1:iter,Energy)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
U = B*Xi;                          % TV pce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
plot(y,F,'g-.')
hold on
plot(y,U_0,'r')
plot(y,U,'b')
xlabel('z')
ylabel('u(0.5,0.5,z)')
% legend('exact','pce','TV-pce','Location','northwest')
legend('exact','pce','TV-pce')
box off

figure
plot(y,abs(F'-U_0),'r')
hold on
plot(y,abs(F'-U),'b')
xlabel('z')
ylabel('Absolute error')
set(gca, 'YScale', 'log')
legend('standard pce','TV pce','Location','northwest')
box off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean & variance
% [mean_mc,var_mc] = Monte_Carlo(f,1e7);
[mean_ex, var_ex] = Exact_Moments(f,l,u,N);
mean_tv = Xi(1);
var_tv = H*(Xi(2:end).^2);

disp('%%%%%%%%%%%%%%%%%%%%%')
disp('Mean')
fprintf('Exact mean: %0.4f \n', mean_ex)
fprintf('Standard PCE Error: %0.4e \n', abs(mean_ex-mean_pce))
fprintf('TV - PCE Error: %0.4e \n', abs(mean_ex-mean_tv))
disp('%%%%%%%%%%%%%%%%%%%%%')
disp('Var')
fprintf('Exact Variance: %0.4f \n', var_ex)
fprintf('Standard PCE Error: %0.4e \n', abs(var_ex-var_pce))
fprintf('TV - PCE Error: %0.4e \n', abs(var_ex-var_tv))
disp('%%%%%%%%%%%%%%%%%%%%%')
disp('L-2 Error')
fprintf('Standard PCE: %0.4e \n', norm(U-U_0))
fprintf('TV - PCE: %0.4e \n', norm(U-F'))