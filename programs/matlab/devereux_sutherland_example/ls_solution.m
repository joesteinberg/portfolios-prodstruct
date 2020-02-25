function varargout = ls_solution(varargin);

% Solves for a first-order or second-order accurate state-space solution. 
% The second-order solution is derived using the Lombardo-Sutherland
% method.
%
% For the first-order solution use
% [F1, F2, F3, P1, P2, P3] = ls_solution(A1, A2, A3, BB, NN, nx, ns, nc);

% For the second-order solution use
% [F1, F2, F3, F4, F5, P1, P2, P3, P4, P5] = ls_solution(A1, A2, A3, A4, A5, BB, NN, SS, nx, ns, nc);

% A1-A5 and NN are the coefficient matrices of the approximated model
% BB is the B vector necessary for the Devereux-Sutherland calculation
% SS is the variance matrix of the innovations
% nx is the number of variables in x
% ns is the number of variables in s
% nc is the number of variables in c
%
% F1-F5 and P1-P5 are the state-space coefficient matrices.
%
% Calls h2v(), v2h(), vech(), perm()
% Calls qzdiv() (by Christopher Sims)
%
% Some sections of code relating to the QZ decomposition are taken
% from solab.m written by Paul Klein.
%
% Written by Alan Sutherland, 20 February 2007
%

if nargout == 6
    A1=varargin{1};   
    A2=varargin{2};
    A3=varargin{3};
    BB=varargin{4};
    NN=varargin{5};
    nx=varargin{6};
    ns=varargin{7};
    nc=varargin{8};    
end

if nargout == 10
    A1=varargin{1};   
    A2=varargin{2};
    A3=varargin{3};
    A4=varargin{4};
    A5=varargin{5};
    BB=varargin{6};
    NN=varargin{7};
    SS=varargin{8};
    nx=varargin{9};
    ns=varargin{10};
    nc=varargin{11};    
end    

nz=ns+nx;
nd=nc+nz;
nr=ns+nc;
itr1=2000;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stake=1.000001;

[s,t,qk,z] = qz(A1,A2);            
[s,t,qk,z] = qzdiv(stake,s,t,qk,z); 

a=qk*A1*z;
b=qk*A2*z;

z11=z(1:ns,1:ns);
z12=z(1:ns,ns+1:nr);
z21=z(ns+1:nr,1:ns);
z22=z(ns+1:nr,ns+1:nr);

a11=a(1:ns,1:ns);
a12=a(1:ns,ns+1:nr);
a22=a(ns+1:nr,ns+1:nr);

b11=b(1:ns,1:ns);
b12=b(1:ns,ns+1:nr);
b22=b(ns+1:nr,ns+1:nr);

ti=inv(b22);
ts=ti*a22;

if rank(z11)<ns;
rank(z11)
 error('Invertibility condition violated')
end

if abs(t(ns,ns))>stake*abs(s(ns,ns)) | abs(t(ns+1,ns+1))<stake*abs(s(ns+1,ns+1));
 warning('Wrong number of stable eigenvalues.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CX=qk*A3;
BX=qk*BB;

C1=CX(1:ns,1:nx);
C2=CX(ns+1:nr,1:nx);

B1=BX(1:ns,1);
B2=BX(ns+1:nr,1);

KK=ti*C2;
M1=KK;

for ii=1:itr1
    KK=ts*KK*NN;
    M1=M1+KK;
end    

tst=max(max(abs(KK)));

if tst>0.00000001;
 warning('M1 not converged.');
end

BM1=ti*B2;

P1=-inv(z22')*M1;
P2=-inv(z22')*z12';
P3=-inv(z22')*BM1;

R1=(a11*z21'+a12*z22')*P1;
R2=(a11*z21'+a12*z22')*P2+(a11*z11'+a12*z12');
R3=(a11*z21'+a12*z22')*P3;

D1=(b11*z21'+b12*z22')*P1+C1;
D2=(b11*z21'+b12*z22')*P2+(b11*z11'+b12*z12');
D3=(b11*z21'+b12*z22')*P3+B1;

F1=inv(R2)*(D1-R1*NN);
F2=inv(R2)*D2;
F3=inv(R2)*D3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 6
    varargout{1}=real(F1);   
    varargout{2}=real(F2);   
    varargout{3}=real(F3);   
    varargout{4}=real(P1);   
    varargout{5}=real(P2);   
    varargout{6}=real(P3);   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 10

    ndd=nd*(nd+1)/2;
    nzz=nz*(nz+1)/2;
    itr2=2000;
   
    H=[NN, zeros(nx,ns); F1, F2];
    C=[eye(nx,nx); zeros(ns,nx)];
    W=[eye(nx,nx),zeros(nx,ns);zeros(ns,nx),eye(ns,ns);P1,P2];

    Lh=h2v(nz*nz,nzz);
    Lc=v2h(nzz,nz*nz);
    FH=Lc*kron(H,H)*Lh;

    Lh=h2v(nx*nx,(nx*(nx+1))/2);
    MH=Lc*kron(C,C)*Lh;

    HKC=kron(H,C);
    CKH=kron(C,H);
    PR=perm(nz,nx);
    PSI=Lc*(CKH+HKC*PR');

    Lh=h2v(nz*nz,nzz);
    Lc=v2h(ndd,nd*nd);
    R=Lc*kron(W,W)*Lh;

    SIGMA=vech(SS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    GG=A4*R+A5*R*FH;
    HH=A5*R*MH;

    GX=qk*GG;
    HX=qk*HH;

    G1=GX(1:ns,:);
    H1=HX(1:ns,:);
    G2=GX(ns+1:nr,:);
    H2=HX(ns+1:nr,:);

    KK=ti*G2;
    M2=KK;

    for ii=1:itr2
        KK=ts*KK*FH;
        M2=M2+KK;
    end    

    tst=max(max(abs(KK)));

    if tst>0.00000001;
     warning('M2 not converged.');
     tst
     KK
    end

    P4=-inv(z22')*M2;

    MX=inv(eye(nc)-ts);
    TX1=MX*ts*M2*MH;
    TX2=MX*ti*H2;
    P5=-inv(z22')*(TX1+TX2);

    R4=(a11*z21'+a12*z22')*P4;
    R5=(a11*z21'+a12*z22')*P5;

    D4=(b11*z21'+b12*z22')*P4+G1;
    D5=(b11*z21'+b12*z22')*P5+H1;

    F4=inv(R2)*(D4-R4*FH);
    F5=inv(R2)*(D5-R5-R4*MH);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    varargout{1}=real(F1);   
    varargout{2}=real(F2);   
    varargout{3}=real(F3);   
    varargout{4}=real(F4);   
    varargout{5}=real(F5);   
    varargout{6}=real(P1);   
    varargout{7}=real(P2);   
    varargout{8}=real(P3);
    varargout{9}=real(P4);   
    varargout{10}=real(P5); 
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
