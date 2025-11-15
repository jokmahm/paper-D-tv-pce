function y = Mc_2D(a,t,z1,z2,mod)
switch mod
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        n = length(z1);
        y = zeros(1,n);
        for i = 1:n
            if z2<=0.5 % && z2(i)>=0.4
                y(i) = (z1(i)/10)*exp(t)*(1-a)^(2/3)*exp(-a);
             else
                y(i) = z1(i)*exp(t)*(1-a)^(2/3)*exp(-a);
             end
        end
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
        n = length(z2);
        y = zeros(1,n);
        for i = 1:n
            if z2(i)<=0.5 % && z2(i)>=0.4
                y(i) = (z1/10)*exp(t)*(1-a)^(2/3)*exp(-a);
             else
                y(i) = z1*exp(t)*(1-a)^(2/3)*exp(-a);
             end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        [m,n] = size(z1);
        y = zeros(m,n);
        for i = 1:m
            for j = 1:n
                if z2(i,j)<=0.5 % && z2(i,j)>=0.4
                    y(i,j) = (z1(i,j)/10)*exp(t)*(1-a)^(2/3)*exp(-a);
                else
                    y(i,j) = z1(i,j)*exp(t)*(1-a)^(2/3)*exp(-a);
                end
            end
        end
end