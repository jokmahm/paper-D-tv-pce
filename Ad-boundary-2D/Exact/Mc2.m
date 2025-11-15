function y = Mc2(a,t,z1,z2,mod)
switch mod
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        if z2<=0.5+0.4*sin(2*pi*z1)
            y = (z1/2)*exp(t)*(1-a)^(2/3)*exp(-a);
        else
            y = z1*exp(t)*(1-a)^(4/3)*exp(-a);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
        [m,n] = size(z1);
        y = zeros(m,n);
        for i = 1:m
            for j = 1:n
                if z2(i,j)<=0.5+0.3*sin(5*pi*z1(i,j)) % 0.5+0.4*sin(2*pi*z1(i,j))
                    y(i,j) = (z1(i,j)/2)*exp(t)*(1-a)^(2/3)*exp(-a);
                else
                    y(i,j) = z1(i,j)*exp(t)*(1-a)^(2/3)*exp(-a);
                end
            end
        end
end