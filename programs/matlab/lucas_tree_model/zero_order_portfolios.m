% FIRST-ORDER SOLUTION AND STEADY STATE PORTFOLIO SHARE CALCULATION
% -------------------------------------------------------------------------------------------

% add first-order coefficients for price dummies, excess return and
% consumption differentials
for i=1:NCOUNTRIES
    s = num2str(i);
    ixf = indsF.(['qm_',s]);
    fy(ixf,indsY.(['q_',s])) = 1;
    fxp(ixf,indsX.(['qm_',s])) = -1;
end

for i=1:(NCOUNTRIES-1)
    s = num2str(i);
    ixf = indsF.(['Rx_',s,num2str(NCOUNTRIES)]);
    fy(ixf,indsY.(['R_',s])) = 1;
    fy(ixf,indsY.(['R_',num2str(NCOUNTRIES)])) = -1;
    fy(ixf,indsY.(['Rx_',s,num2str(NCOUNTRIES)])) = -1;
end

for i=1:(NCOUNTRIES-1)
    s = num2str(i);
    ixf = indsF.(['cD_',s,num2str(NCOUNTRIES)]);
    fy(ixf,indsY.(['c_',s])) = gama;
    fy(ixf,indsY.(['c_',num2str(NCOUNTRIES)])) = -gama;
    fy(ixf,indsY.(['cD_',s,num2str(NCOUNTRIES)])) = -1;
end

% evaluate first-order expansion and construct state-space matrices
nfz = zeros(size(fz));
nfx = zeros(size(fx));
nfy = zeros(size(fy));
nfzp = zeros(size(fzp));
nfxp = zeros(size(fxp));
nfyp = zeros(size(fyp));

nfz(:) = eval(fz(:));
nfx(:) = eval(fx(:));
nfy(:) = eval(fy(:));
nfzp(:) = eval(fzp(:));
nfxp(:) = eval(fxp(:));
nfyp(:) = eval(fyp(:));

A1 = -[nfxp, nfyp];
A2 = [nfx, nfy];
A3 = nfz;

% vector of coefficients on the iid variables xi
B=zeros(nx+ny,NCOUNTRIES-1);
for i=1:(NCOUNTRIES-1)
    s = num2str(i);
    B(indsF.(['bc_',s]),i)=1;
end

% first-order solution using Lombardo-Sutherland algorithm
[F1, F2, F3, P1, P2, P3] = ls_solution_js(A1, A2, A3, B, N, nz, nx, ny);

% covert to Devereux-Sutherland format
[Q, M, nq] = ls2ds_js(F1, F2, P1, P2, P3, N, 1);

% R1, R2, D1, D2 matrices
R1=[];
R2=[];
D1=[];
D2=[];
for i=1:(NCOUNTRIES-1)
    s=num2str(i);
    R1 = [R1;Q{indsY.(['Rx_',s,num2str(NCOUNTRIES)]),1}];
    R2 = [R2;Q{indsY.(['Rx_',s,num2str(NCOUNTRIES)]),2}];
    D1 = [D1;Q{indsY.(['cD_',s,num2str(NCOUNTRIES)]),1}];
    D2 = [D2;Q{indsY.(['cD_',s,num2str(NCOUNTRIES)]),2}];
end

% symmetric initial guess
alpha0 = exp(q_1)/NCOUNTRIES * ones(NCOUNTRIES-1);
for i=1:(NCOUNTRIES-1)
    alpha0(i,i) = alpha0(i,i) - exp(q_1);
end
%alpha0 = zeros(NCOUNTRIES-1);
%for i=1:(NCOUNTRIES-1)
%    for j=1:(NCOUNTRIES-1)
%        if(i~=j)
%            alpha0(i,j) = exp(q_1)/(NCOUNTRIES*theta);
%        end
%    end
%end
%for i=1:(NCOUNTRIES-1)
%    alpha0(i,i) = -sum(alpha0(:,i))-exp(q_1);
%end

% solve the E_t[cD*Rx] equation
params = {};
params.NCOUNTRIES = NCOUNTRIES;
params.R1 = R1;
params.R2 = R2;
params.D1 = D1;
params.D2 = D2;
params.SIGMA = SIGMA;

options = optimset('fsolve');
options = optimset('TolFun',1e-10);
options = optimset('TolX',1e-10);
options = optimset('Display','off');
[alpha,fval,exitflag] = fsolve( @(alpha) share_func(alpha,params) , alpha0, options);
if(exitflag ~= 1)
    warning('Warning! Fsolve failed to find equilibrium portfolio shares!');
end

% calculate alpha and lambda
alpha2 = [alpha*beta,zeros(NCOUNTRIES-1,1)];
for i=1:(NCOUNTRIES-1)
    s=num2str(s);
    alpha2(i,NCOUNTRIES) = eval(['WEALTH_',s]) - sum(alpha2(i,:));
end

lambda = zeros(size(alpha2));
for i=1:(NCOUNTRIES-1)
    for j=1:NCOUNTRIES
        tmp = alpha2(i,j);
        if(i==j)
            tmp = tmp + exp(eval(['q_',num2str(i)]));
        end
        lambda(i,j) = tmp/exp(eval(['q_',num2str(j)]));
    end
end

% print home and foreign shares
fprintf('\nModel home share: %f\n',lambda(1,1));
fprintf('Model foreign shares: ');
fprintf('%f ',lambda(1,2:(NCOUNTRIES)));
fprintf('\n');

lambda_th = ones(NCOUNTRIES-1,NCOUNTRIES) / (theta*NCOUNTRIES);
for i=1:(NCOUNTRIES-1)
    lambda_th(i,i) = (1-NCOUNTRIES+NCOUNTRIES*theta)/(NCOUNTRIES*theta);
end
fprintf('\nTheoretical home share: %f\n',lambda_th(1,1));
fprintf('Theoretical foreign shares: ');
fprintf('%f ',lambda_th(1,2:(NCOUNTRIES)));
fprintf('\n');