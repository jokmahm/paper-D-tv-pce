function y = Mc7_2D(a,t,z1,z2,mod)
switch mod
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        if z2<=0.5+0.1*sin(2*pi*z1)
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
                if z2(i,j)<=0.5 && z1(i,j)>=0 && z1(i,j)<0.35
                    
                    y(i,j) = (0.5)*exp(t)*(1-a)^(2/3)*exp(-a);      
                    
                elseif z2(i,j)<=0.5 && z1(i,j)>=0.35 && z1(i,j)<0.65
                    
                    y(i,j) = (1)*exp(t)*(1-a)^(2/3)*exp(-a);
                    
                elseif z2(i,j)<=0.5 &&  z1(i,j)>=0.65 && z1(i,j)<=1
                    
                    y(i,j) = (1.5)*exp(t)*(1-a)^(2/3)*exp(-a);
                    
                elseif z2(i,j)>0.5 && z1(i,j)>=0 && z1(i,j)<0.35
                    
                    y(i,j) = 2*(0.5)*exp(t)*(1-a)^(2/3)*exp(-a);      
                    
                elseif z2(i,j)>0.5 && z1(i,j)>=0.35 && z1(i,j)<0.65
                    
                    y(i,j) = 2*(1)*exp(t)*(1-a)^(2/3)*exp(-a);
                    
                elseif z2(i,j)>0.5 &&  z1(i,j)>=0.65 && z1(i,j)<=1
                    
                    y(i,j) = 2*(1.5)*exp(t)*(1-a)^(2/3)*exp(-a);
                    
                end
            end
        end
end