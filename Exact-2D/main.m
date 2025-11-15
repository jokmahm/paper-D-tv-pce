clc; clear;
% Recurrance form of Legender construction
% Exact derivative of the Legender polynomials
l = 0;                           % parameter uncertainty lower bound
u = 1;                           % parameter uncertainty upper bound
mod = 3;

global idx
% profile on
lambda = 250;
p = 30;                          % degree of the polynomial

n = 20;                          % number of discritization
N = p+10;                        % order Gauss quaderature integration rule
x1 = linspace(l,u,n);            % discritization for simpson integration
x2 = linspace(l,u,n);
[X1, X2] = meshgrid(x1,x2);

% f = @(x,y) exp(-abs(x-y));         
% f = @(x,y) sin(x+y);
% f = @(x,y) Genze6(0.5,x,y,4);
% f = @(x,y) Test_2D(x,y);

f = @(x,y) Mc_2D(0.5,0.5,x,y,mod);   % Exact solution
% f = @(x,y) Mc1_2D(0.5,0.5,x,y);

%g = @(x,y) Upwind(x,y,mod);          % numerical solution of the PDE
% Y = zeros(n,n);
% for i= 1:n
%     for j=1:n
%         Y(i,j) = g(X1(i,j),X2(i,j));
%     end
% end

Z = f(X1,X2);
Q = zeros(size(Z));
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

B = basis(p,l,u,X1,X2);
[B_p, B_x, B_y, W, V] = Guass_quad_integration(p,l,u,l,u,N,f);

for i = 1:Np
   I = B_p(:,:,i).*V.*W;
   b(i) = (u-l)*(u-l)/4*sum(sum(I));
end 

gamma = Mass(p,l,u);              % orthogonality coefficients
M = diag(gamma);                  % Mass matrix
Xi = M\b;                         % coefficients
c = Xi;
r = lambda*b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard mean and variance report
mean_pce = Xi(1);
C = Xi(2:end);
H = gamma(2:end);
var_pce = (C.^2)'*H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = zeros(size(Z));
for i = 1:Np
    G = G+Xi(i)*B(:,:,i);          % solution of the standard PCE
end

iter = 20;
Energy = zeros(1,iter);
for k=1:iter                       % non-linear loop
    fprintf('Iteration %d \n', k)         % f', k, norm(r))
    K = Assembler2D(p, Xi, B_x, B_y, W, l, u, l, u);
    J = K+lambda*M;
    Xi = J\r;                      % solve for correction
    Q = zeros(n,n);
    for i = 1:Np
        Q = Q+Xi(i)*B(:,:,i);
    end
    [ux,uy] = gradu(Xi,p,B_x,B_y);
    Energy(k) = norm(sqrt(ux.^2+uy.^2),1)+0.5*lambda*norm(Q-Z,2);
    fprintf('Energy Norm is %0.4f \n', Energy(k));
end

figure
plot(1:iter,Energy)

%%%%
n = 19;                          % number of discritization
x1 = linspace(l,u,n);             % discritization for simpson integration
x2 = linspace(l,u,n);
[X1, X2] = meshgrid(x1,x2);
Z = f(X1,X2);
B = basis(p,l,u,X1,X2);
G = zeros(n,n);
for i = 1:Np
    G = G+c(i)*B(:,:,i);          % solution of the standard PCE
end
Q = zeros(n,n);
for i = 1:Np
    Q = Q+Xi(i)*B(:,:,i);
end
Er1 = abs(Z-G);
Er2 = abs(Z-Q);
%%%%

figure
subplot(2,2,1), surf(X1,X2,Z)
xlabel('x')
ylabel('y')
zlabel('Exact')
subplot(2,2,2), surf(X1,X2,G)
xlabel('x')
ylabel('y')
zlabel('Standard PCE')
subplot(2,2,3), surf(X1,X2,Q)
xlabel('x')
ylabel('y')
zlabel('TV result')
subplot(2,2,4), surf(X1,X2,Er2)
xlabel('x')
ylabel('y')
zlabel('Error')
set(gca, 'ZScale', 'log')

figure()
subplot(1,2,1), surf(X1,X2,Er1)
xlabel('z1')
ylabel('z2')
zlabel('PCE error')
set(gca, 'ZScale', 'log')
subplot(1,2,2), surf(X1,X2,Er2)
xlabel('z1')
ylabel('z2')
zlabel('TV-PCE error')
set(gca, 'ZScale', 'log')
 
figure
plot(x1,Er1(:,10))
hold on
plot(x1,Er2(:,10))
set(gca, 'YScale', 'log')
xlabel('z1')
ylabel('absolute error')
legend('Standard-pce','TV-pce')
box off

figure
plot(x1,sqrt(sum((Z-G).^2,2)))
hold on
plot(x1,sqrt(sum((Z-Q).^2,2)))
set(gca, 'YScale', 'log')
xlabel('z1')
ylabel('L-2 error')
legend('Standard-pce','TV-pce')
box off

figure
plot(x2,sqrt(sum((Z-G).^2,1)))
hold on
plot(x2,sqrt(sum((Z-Q).^2,1)))
set(gca, 'YScale', 'log')
xlabel('z2')
ylabel('L-2 error')
legend('Standard-pce','TV-pce')
box off

% figure
% subplot(1,3,1), surf(X1,X2,Z)
% xlabel('x')
% ylabel('y')
% zlabel('Exact')
% subplot(1,3,2), surf(X1,X2,G)
% xlabel('x')
% ylabel('y')
% zlabel('Standard PCE')
% subplot(1,3,3), surf(X1,X2,Q)
% xlabel('x')
% ylabel('y')
% zlabel('Total Variation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mean & variance
% [mean_mc,var_mc] = Monte_Carlo(f,1e7);
[mean_ex, var_ex] = Exact_Moments(f,l,u,l,u,N);
mean_tv = Xi(1);
C = Xi(2:end);
var_tv = (C.^2)'*H;

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
fprintf('Standard PCE: %0.4e \n', norm(Z-G))
fprintf('TV - PCE: %0.4e \n', norm(Z-Q))

% profile viewer
% profile off