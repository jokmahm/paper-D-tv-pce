clc; clear;
% Recurrance form of Legender construction
% Exact derivative of the Legender polynomials
% Add Numann boundary condition
global idx
% profile on
p = 31;                             % degree of the polynomial
l1 = 0.5;                           % parameter uncertainty lower bound
u1 = 1.5;                           % parameter uncertainty upper bound
l2 = 0;
u2 = 1;
lambda = 10*p;

n = 20;                           % number of discritization
N = p+10;                         % order Gauss quaderature integration rule
x1 = linspace(l1,u1,n);           % discritization for simpson integration
x2 = linspace(l2,u2,n);
[X1, X2] = meshgrid(x1,x2);

% f = @(x,y) sin(x+y);
%f = @(x,y) Mc_2D(0.5,0.5,x,y,3);
f = @(x,y) Mc2(0.5,0.5,x,y,2);

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

B = basis(p,l1,u1,l2,u2,X1,X2);
[B_p, B_x, B_y, W, V] = Guass_quad_integration(p,l1,u1,l2,u2,N,f);

for i = 1:Np
   I = B_p(:,:,i).*V.*W;
   b(i) = ((u1-l1)/2)*((u2-l2)/2)*sum(sum(I));
end 

gamma = Mass(p,l1,u1,l2,u2);              % orthogonality coefficients
M = diag(gamma);                  % Mass matrix
Xi = M\b;                         % coefficients
c = Xi;

G = zeros(size(Z));
for i = 1:Np
    G = G+c(i)*B(:,:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add line integral
L = zeros(Np,1);                % boundary values
[Bx0, Bx1, B0y, B1y, w] = Linex(p,l1,u1,l2,u2,N);
%h1 = f(1.5*ones(1,N),linspace(0,1,N))./sqrt(f(1.5*ones(1,N),linspace(0,1,N)).^2+1e-4);
%h2 = f(0.5*ones(1,N),linspace(0,1,N))./sqrt(f(0.5*ones(1,N),linspace(0,1,N)).^2+1e-4);

for i = 1:Np 
       % line integral from exact function 
%    I1 = -((fy_x0./gx0).*Bx0(i,:).*w);
%    I2 = ((fx_1y./g1y).*B1y(i,:).*w);
%    I3 = -((fy_x1./gx1).*Bx1(i,:).*w);
%    I4 = ((fx_0y./g0y).*B0y(i,:).*w);

   I1 = 0;
   I2 = -(0.7.*B1y(i,:).*w); % -(0.33.*B1y(i,:).*w);
   I3 = 0;
   I4 = (0.8.*B0y(i,:).*w); %(0.46.*B0y(i,:).*w);
   
   L(i) = ((u1-l1)/2)*(sum(I1)+sum(I3)+sum(I2)+sum(I4));
end
r = lambda*b-L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r = lambda*b;

iter = 50;
Energy = zeros(1,iter);
for k=1:iter                     % non-linear loop
    fprintf('Iteration %d \n', k)         % f', k, norm(r))
    K = Assembler2D(p, Xi, B_x, B_y, W, l1, u1, l2, u2);
    J = K+lambda*M;
    
    Xi = J\r;                    % solve for correction
    Q = zeros(n,n);
    for i = 1:Np
        Q = Q+Xi(i)*B(:,:,i);
    end
    [u_x, u_y] = gradu(Xi,p,B_x,B_y);
    Energy(k) = norm(sqrt(u_x.^2+u_y.^2),1)+0.5*lambda*norm((Q-Z),2)^2;
    fprintf('Energy Norm is %0.4f \n', Energy(k));
end
figure
plot(1:iter,Energy)

%%%%
n = 50;                          % number of discritization
x1 = linspace(l1,u1,n);             % discritization for simpson integration
x2 = linspace(l2,u2,n);
[X1, X2] = meshgrid(x1,x2);
Z = f(X1,X2);
B = basis(p,l1,u1,l2,u2,X1,X2);
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
