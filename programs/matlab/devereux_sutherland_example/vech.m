function X = vech(Y);

[ng1,ng2]=size(Y);

X=zeros(ng2*(ng2+1)/2,1);

kk=1;

for jj=1:ng2;
    for ii=1:jj; 
        X(kk,1)=Y(ii,jj);
        kk=kk+1;
    end
end    
