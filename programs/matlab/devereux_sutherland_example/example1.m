% solves example model

% calls modeval(), ls_solution(), ls2ds()

clear all

load 'example.mat' model;

% model dimensions

nx=4;
ns=3;
nc=10;

% parameter values

beta=0.98;
rho=2.0;
zeta_y=0.9;
zeta_m=1.0;

[A1,A2,A3,A4,A5,NN,SIGMA]=modeval(model,'BT',beta,'RH',rho,'DY',zeta_y,'DM',zeta_m);

SIGMA=[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

% row indices for W, rh, rf, rx, Cd

iW=1;

irh=5;
irf=6;

irx=9;
icd=10;

% set up B vector

BB=zeros(ns+nc,1);
BB(iW,1)=1;

% first-order solution using Lombardo-Sutherland algorithm

[F1, F2, F3, P1, P2, P3] = ls_solution(A1, A2, A3, BB, NN, nx, ns, nc);

% covert to Devereux-Sutherland format

[Q, M, nq] = ls2ds(F1, F2, P1, P2, P3, NN, 1);

% R1, R2, D1, D2 matrices

R1=Q{irx,1};
R2=Q{irx,2};

D1=Q{icd,1};
D2=Q{icd,2};

% calculate alpha-tilde

alpha_tilde=inv(R2*SIGMA*D2'*R1'-R2*SIGMA*R2'*D1)*R2*SIGMA*D2';

% column indices to identify where alpha-tilde appears in the A-tilde matrices

irhrh=(nx+ns+irh)*((nx+ns+irh)+1)/2;
irfrf=(nx+ns+irf)*((nx+ns+irf)+1)/2;
irfW=(nx+ns+irf-1)*(nx+ns+irf)/2+nx+iW;

% form A-tilde matrices

At1=A1;
At2=A2;
At3=A3;
At4=A4;
At5=A5;

At2(iW,irh+ns)=alpha_tilde;  
At2(iW,irf+ns)=-alpha_tilde;   

At4(iW,irhrh)=alpha_tilde/2;
At4(iW,irfrf)=-alpha_tilde/2;
At4(iW,irfW)=1/beta;

% second-order solution using Lombardo-Sutherland algorithm

[Ft1, Ft2, Ft3, Ft4, Ft5, Pt1, Pt2, Pt3, Pt4, Pt5] = ls_solution(At1, At2, At3, At4, At5, BB, NN, SIGMA, nx, ns, nc);

% covert to Devereux-Sutherland format

[Q, M, nq] = ls2ds(Ft1, Ft2, Pt1, Pt2, Pt3, Pt4, NN, 1);

% D-tilde and R-tilde matrices

Dt1=Q{icd,1};
Dt2=Q{icd,2};
Dt5=Q{icd,5};

Rt2=Q{irx,2};
Rt5=Q{irx,5};

% calculate gamma

gamma_prime=-inv(Rt2*SIGMA*Rt2'*Dt1)*(Rt2*SIGMA*Dt5'+Dt2*SIGMA*Rt5');
gamma=gamma_prime';

% end

