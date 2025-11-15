function J = Assembeler(Xi,B_x,W,p,a,b)
du = abs(B_x*Xi); 
du(du==0) = 1e-8;
J = zeros(p+1,p+1);
for i = 1:p+1
    for j = 1:p+1
        h = (1./du).*B_x(:,i).*B_x(:,j);
        E = h.*W';
        J(i,j) = (b-a)/2*sum(E);   % gauss integration rule
    end
end