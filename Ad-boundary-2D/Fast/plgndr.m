function pn=plgndr(z,n)

pn=zeros(length(z),n+1);

% code if(n==0)

if (n==0)

pn(:,1)=1;
return
end

pnm = flgndr(z,n);
pn(:,n+1)=pnm(:,n+1);

pnm = flgndr(z,n-1);
pn(:,n)=pnm(:,n);

for nc=n-1:-1:1

pn(:,nc)=((2*nc+1)*z.*pn(:,nc+1)- (nc+1)*pn(:,nc+2))/nc ;
end