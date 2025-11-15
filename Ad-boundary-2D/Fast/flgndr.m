function pnm = flgndr(z,n)
nz = length(z);
pnm = zeros(nz,n+1);

fac = prod(2:n);
sqz2 = sqrt((1.0-z.*z));
sqz2(sqz2==0) = 1e-1;
hsqz2 = 0.5*sqz2;
ihsqz2 = z./hsqz2;

if(n==0)

pnm(:,1) = 1.0;
return
end
if(n==1)
pnm(:,1) = -0.5*sqz2;
pnm(:,2) = z;
pnm(:,3) = sqz2;
return
end
pnm(:,1) = (1-2*abs(n-2*floor(n/2)))*hsqz2.^n/fac;
pnm(:,2) = -pnm(:,1)*n.*ihsqz2;
for mr=1:2*n-1

pnm(:,mr+2)=(mr-n).*ihsqz2.*pnm(:,mr+1)-(2*n-mr+1)*mr*pnm(:,mr);
end
end