function y = Mc_2D(a,t,z1,z2,mod)
[m,n] = size(z1);
y = zeros(m,n);
u = exp(t)*(1-a)^(2/3)*exp(-a);
switch mod
    case 1       
        for i = 1:m
            for j = 1:n
                if z2(i,j)<=0.5 
                    y(i,j) = (1/2)*u;
                else
                    y(i,j) = u;
                end
            end
        end
    case 2
        for i = 1:m
            for j = 1:n
                if z2(i,j)<=0.5+0.1*sin(2*pi*z1(i,j))
                    y(i,j) = (1/2)*u; 
                else
                    y(i,j) = u; 
                end
            end
        end
    case 3
        r1 = 0.8;
        r2 = 0.8;
        for i = 1:m
            for j = 1:n
                if (z1(i,j)-r1)^2+(z2(i,j)-r2)^2 <= 1/10
                    y(i,j) = u;
                else
                    y(i,j) = (1/2)*u;
                end
            end
        end
end