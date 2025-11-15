function y = Genze6(x,alpha,beta,mod)
switch mod
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        if beta<=x 
            y = x-1;
        else
            y = exp(alpha*x);
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2
        nx = length(x);
        y = zeros(nx,1);
        for i = 1:nx
            if beta<=x(i) 
                y(i) = x(i)-1;
            else
                y(i) = exp(alpha*x(i));
            end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 3
        na = length(alpha);
        nb = length(beta);
        y = zeros(na,nb);
        for i = 1:na
            for j = 1:nb
                if beta(j)<=x 
                    y(i,j) = x-1;
                else 
                    y(i,j) = exp(alpha(i)*x);
                end
            end
        end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 4
        [na, nb] = size(alpha);
        y = zeros(na,nb);
        for i = 1:na
            for j = 1:nb
                if beta(i,j)<=x 
                    y(i,j) = x-1;
                else 
                    y(i,j) = exp(alpha(i,j)*x);
                end
            end
        end
end