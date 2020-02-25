%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% finish setting up numerical version of linearized system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    fy(ixf,indsY.(['fc_',s])) = -gama;
    fy(ixf,indsY.(['pfc_',s])) = -1;
    fy(ixf,indsY.(['fc_',num2str(NCOUNTRIES)])) = gama;
    fy(ixf,indsY.(['pfc_',num2str(NCOUNTRIES)])) = 1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve for decision rules
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create initial guess for portfolio shares
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note that country i's share of total wealth in the economy is
% sharei = (NFA_i + MKTCAP_i)/sum(MKTCAP_j)
% where NFA_i is net foreign assets and MKTCAP_i is the stock
% market capitalization of country i. NFA_j does not show up on the
% x-axis because sum(NFA_j) = 0.
%
% start by assuming

% total world stock market cap = sum(q_j)
wstk = 0;
for i=1:NCOUNTRIES
    s=num2str(i);
    wstk = wstk + exp(eval(['q_',s]));
end

% start by assuming each country i holds the same i-specific sharei
% of each country j's stock, where sharei equals
% country i's stock market cap + NFA as a fraction of world stock market cap
% then change country i's holdings of its own stock so that the sum
% of its alphas equals WEALTH_i
if(sens==0 && (outer==1 || outer==7))
    ALPHA0 = zeros(NCOUNTRIES);
    for i=1:NCOUNTRIES
        s=num2str(i);
        sharei = (exp(eval(['q_',s]))+eval(['WEALTH_',s]))/wstk;
        for j=1:NCOUNTRIES
            s2=num2str(j);
            ALPHA0(i,j) = sharei * exp(eval(['q_',s2]));
        end
        ALPHA0(i,i) = ALPHA0(i,i) - exp(eval(['q_',s]));
    end

    % check that everything adds up
    for i=1:NCOUNTRIES
        s=num2str(i);
        tmp = eval(['WEALTH_',s]) - sum(ALPHA0(i,:));
        if abs(tmp)>TINY
            error(['alpha0(',s,':) do not sum to WEALTH(',s,')!']);
        end
        tmp = sum(ALPHA0(:,i));
        if abs(tmp)>TINY
            error(['alpha0(',s,':) do not sum to zero!']);
        end
    end
    % take the (I-1) x (I-1) sub-matrix that we actually need to solve for
    ALPHA0 = ALPHA0(1:(NCOUNTRIES-1),1:(NCOUNTRIES-1));
else
    ALPHA0 = ALPHA;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solve for steady-state portfolios (and wedges if necessary)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

params = {};
params.NCOUNTRIES = NCOUNTRIES;
params.R1 = R1;
params.R2 = R2;
params.D1 = D1;
params.D2 = D2;
params.SIGMA = SIGMA;
params.TINY = TINY;
params.DIV = DIV;
params.ALPHA0 = ALPHA0;
params.beta = beta;
params.WEALTH_1 = WEALTH_1;
params.q_1 = q_1;
params.WEALTH_2 = WEALTH_2;
params.q_2 = q_2;
params.WEALTH_3 = WEALTH_3;
params.q_3 = q_3;
params.WEALTH_4 = WEALTH_4;
params.q_4 = q_4;
params.SIGMA = SIGMA;

% if we are in the 1995 benchmark, we need to first solve for the wedges
% to match the diversification data
% solve E_t[tau*cD*Rx] equation, taking wedges as given, and
% retrieve data
if(tauflag==1)

    if(sens==0)

        % first solve for portfolios with no wedges so we don't start too far away
        [F,ALPHA,ALPHA3,lambda] = zero_order_portfolios2(TAU,params);

        options = optimset('fsolve');
        options = optimset(options,'TolFun',1e-2);
        options = optimset(options,'TolX',0.0);
        options = optimset(options,'Display','iter');
        options = optimset(options,'TypicalX',0.01*diag(SIGMA));
        options = optimset(options,'Algorithm','trust-region-dogleg');
        [TAU,fval,exitflag] = fsolve( @(TAU) zero_order_portfolios2(TAU,params) , TAU, options);
        if(max(abs(fval)) > 1e-1)
            error('Warning! Fsolve failed to find wedges!');
        end

        %if(outer<7)
        %    TAU0=TAU;
        %else
        %    TAU1=TAU;
        %end
    end
end

if(outer<7)
    TAU=TAU0;
else
    TAU=TAU1;
end

% solve E_t[tau*cD*Rx] equation, taking wedges as given, and
% retrieve data
[F,ALPHA,ALPHA3,lambda] = zero_order_portfolios2(TAU,params);

% print home and foreign shares
%if verbose==1
%fprintf('\n\tPortfolio weights:\n');
%for i=1:NCOUNTRIES
%    fprintf('\t%f ',ALPHA3(i,:));
%    fprintf('\n');
%end
%fprintf('\n');
    %end
fprintf('\n\t\tDiversification:');
for i=1:NCOUNTRIES
    fprintf([' ',num2str(100*(1-ALPHA3(i,i)))]);
end
fprintf('\n');

% write them to text file
% print home and foreign shares
if sens==0
    fid = fopen([output_path,'results',str,mname,'.txt'],'wb');
else
    fid = fopen([output_path,'results-s',num2str(sens),str,mname,'.txt'],'wb');
end

for i=1:NCOUNTRIES
    fprintf(fid,'%f ',ALPHA3(i,:));
    fprintf(fid,'\n');
end
fprintf('\n');
fclose(fid);
