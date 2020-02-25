function X = h2v(nr,nc);

X=zeros(nr,nc);

nn=2*nc-nr;

ir=1;
ic=1;

for ii=1:nr;  
    
    if ir<=ic
        jj=(ic*(ic-1)/2)+ir;
    else
        icc=ir;
        irr=ic;
        jj=(icc*(icc-1)/2)+irr;
    end  
    
    X(ii,jj)=1;
    
    ir=ir+1;
    if ir>nn
        ic=ic+1;
        ir=1;
    end 
    
end    
