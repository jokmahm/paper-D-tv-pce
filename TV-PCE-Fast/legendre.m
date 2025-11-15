function l = legendre(n,x)
% HERMITE: compute the legendre polynomials.
% 
%   l = legendre(n)
%   l = legendre(n,x)
% 
% Inputs:
%   - n is the order of the Hermite polynomial (n>=0).
%   - x is (optional) values to be evaluated on the resulting legendre
%     polynomial function.
% 

if( n<0 ), error('The order of legendre polynomial must be greater than or equal to 0.'); end
% again check n is an integer
if( 0~=n-fix(n) ), error('The order of legendre polynomial must be an integer.'); end
% call the legendre recursive function.
l = legendre_rec(n);
% evaluate the legendre polynomial function, given x
if( nargin==2 )
    y = l(end) * ones(size(x));
    p = 1;
    for i=length(l)-1:-1:1
        y = y + l(i) * x.^p;
        p = p+1;
    end
    
    % restore the shape of y, the same as x
    l = reshape(y,size(x));
end
function l = legendre_rec(n)
% This is the reccurence construction of a legendre polynomial, i.e.:
%   L0(x) = 1
%   L1(x) = x
%   L[n+1](x) = [(2n+1)/(n+1)] Ln(x) - [n/(n+1)] L[n-1](x)
if( n==0 ), l = 1;
elseif( n==1 ), l = [1 0];
else
    
    h1 = zeros(1,n+1);
    h1(1:n) = ((2*n-1)/n)*legendre_rec(n-1);
    
    h2 = zeros(1,n+1);
    h2(3:end) = ((n-1)/n)*legendre_rec(n-2);
    
    l = h1 - h2;
    
end