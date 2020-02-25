function X = v2h(nr,nc);

X=zeros(nr,nc);

nn=2*nr-nc;

ir=1;
ic=1;

for ii=1:nr;  
    
    jj=(ic-1)*nn+ir;
    
    X(ii,jj)=1;
    
    ir=ir+1;
    if ir>ic
        ic=ic+1;
        ir=1;
    end 
    
end    
