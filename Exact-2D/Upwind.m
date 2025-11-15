function y = Upwind(z1,z2,mod)
% upwind scheme to compute numerical solution 
% of Mc-kenderik von forster equation
A = 1;
T = 1;

a = 0:0.01:A;
t = 0:0.01:T;

Nt = length(t);           % number of time discretization
Na = length(a);           % number of age discretization

dt = T/Nt;
da = A/Na;

u = zeros(Nt, Na);

u(1,:) = initial_cond(a);              % initial condition
for i = 1:Nt-1
    c = Boundary(a).*u(i,:);
    u(i+1,1) = da*trap(c);    % boundary condition
    for j = 2:Na-1
       u(i+1,j) =  u(i,j)-(dt/da)*(u(i,j)-u(i,j-1))-dt*(mortality(a(j))+Jump(a(j),z1,z2,mod))*u(i,j); 
    end
end    

y = u(51,51);