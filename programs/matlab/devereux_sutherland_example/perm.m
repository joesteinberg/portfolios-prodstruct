function P = perm(n,m);
nx=n*m;
P=zeros(nx,nx);
ii=1;
jj=0;
for j=1:nx
    jj=jj+1;
    if (jj>m)
        jj=1;
        ii=ii+1;
    end 
    i=ii+n*(jj-1);
    P(i,j)=1;
end
