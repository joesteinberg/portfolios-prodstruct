function [Q, M, nq] = ls2ds(varargin);

% Converts output from Lombardo-Sutherland solution routine
% into the format required for Devereux-Sutherland calculations.
%
% For the first-order solution use 
% [Q, M, nq] = ls2ds(F1, F2, P1, P2, P3, NN, ird);
%
% For the second-order solution use
% [Q, M, nq] = ls2ds(F1, F2, P1, P2, P3, P4, NN, ird);
%
% F1, F2, P1, P2, P3, P4, NN are the coefficients from the
% Lombardo-Sutherland solution. 
%
% Q is a (nc X 6) cell array where nc is the number of elements in c.
% Row i of Q gives the coefficients in the solution for variable i in c.
% Only the solutions for the variables in c are returned.
% Q{i,1} is the coefficient on xi 
% Q{i,2} is the (1 X nx) row vector of coefficients on epsilon 
% Q{i,3} is the (1 X nz) row vector of coefficients on z
% Q{i,4} is the (nx X nx) matrix of coefficients on epsilon.epsilon'
% Q{i,5} is the (nz X nx) matrix of coefficients on z.epsilon'
% Q{i,6} is the (nz X nz) matrix of coefficients on z.z'
% where nx is the number of variables in x
% and nz is the number of variables in z.
%
% M is a (nz X nq) matrix of coefficients which map the z vector
% onto the first nq variables in z
%
% ird=1 indicates that M and nq should be calculated
% ird=0 indicates that M and nq should not be calculated 
% in which case M is set to unity and nq is set to nz
%
% Written by Alan Sutherland, 20 February 2007
%

if nargin == 7
    F1=varargin{1};   
    F2=varargin{2};
    P1=varargin{3};
    P2=varargin{4};
    P3=varargin{5};
    NN=varargin{6};
    ird=varargin{7};
end

if nargin == 8
    F1=varargin{1};   
    F2=varargin{2};
    P1=varargin{3};
    P2=varargin{4};
    P3=varargin{5};
    P4=varargin{6};
    NN=varargin{7};
    ird=varargin{8};
end    

[nc, nx]=size(P1);
[nc, ns]=size(P2);

nz=ns+nx;

if ird==1

    XX=[NN zeros(nx,ns); F1*NN F2];
    CC=[eye(nx,nx); F1];

    ev=eig(XX);
    pe=find(ev>0.00000001);
    nq=length(pe);
    nzero=nz-nq;

    pick=[eye(nq,nq) zeros(nq,nzero)];

    [GS,pickn,ok,uis,uiu,vs]=gstate(XX,CC,pick);

    AX=GS(1:nq,:);
    BZ=GS(nq+1:nz,:)*inv(AX);
    
    CC1=CC(1:nq,:);
    CC2=CC(nq+1:nz,:);

    tst=max(max(abs(BZ*CC1-CC2)));
    if tst>0.00000001;
        warning('M invalid.');
    end
    
    M=[eye(nq,nq); BZ];

else    
    
    M=1;
    nq=nz;
    nzero=0;
    
end   

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp3=([P1*NN P2]*M);

for i=1:nc
    PA{i,1}=P3(i,1);
    PA{i,2}=P1(i,:);
    PA{i,3}=tmp3(i,:);
end 

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 8

    nzz=nz*(nz+1)/2;
    nqq=nq*(nq+1)/2;
    nxx=nx*(nx+1)/2;
    nxq=nx*nq;
    
    U1=[NN zeros(nx,ns); zeros(ns,nx) eye(ns,ns)]*M;
    U2=[eye(nx,nx); zeros(ns,nx)];

    Lc=v2h(nzz,nz*nz);
    Lh=h2v(nx*nx,(nx*(nx+1))/2);
    X1=Lc*kron(U2,U2)*Lh;

    U12=kron(U1,U2);
    U21=kron(U2,U1);
    PR=perm(nq,nx);
    X2=Lc*(U21+U12*PR');

    Lh=h2v(nq*nq,nqq);
    X3=Lc*kron(U1,U1)*Lh;

    TT4=P4*X1;
    TT5=P4*X2;
    TT6=P4*X3;

    for i=1:nc
    
        ii=1;
        jj=1;
        tmp4=zeros(nx,nx);
        for k=1:nxx;
            tmp4(ii,jj)=tmp4(ii,jj)+TT4(i,k)/2;
            tmp4(jj,ii)=tmp4(jj,ii)+TT4(i,k)/2;
            ii=ii+1;
            if ii>jj;
                ii=1;
                jj=jj+1;
            end    
        end
        PA{i,4}=tmp4;

        ii=1;
        jj=1;
        for k=1:nxq;
            tmp5(ii,jj)=TT5(i,k);
            ii=ii+1;
            if ii>nq;
                ii=1;
                jj=jj+1;
            end    
        end
        PA{i,5}=tmp5;
    
        ii=1;
        jj=1;
        tmp6=zeros(nq,nq);
        for k=1:nqq;
            tmp6(ii,jj)=tmp6(ii,jj)+TT6(i,k)/2;
            tmp6(jj,ii)=tmp6(jj,ii)+TT6(i,k)/2;
            ii=ii+1;
            if ii>jj;
                ii=1;
                jj=jj+1;
            end    
        end
        PA{i,6}=tmp6;
    
    end 

end

Q=PA;
