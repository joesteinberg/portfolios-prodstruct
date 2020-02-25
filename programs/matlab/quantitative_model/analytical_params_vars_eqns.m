% *****************************************************************************
% declare symbolic variables
% *****************************************************************************

% ---------------------------------------------------------
% parameters
% ---------------------------------------------------------

% gross output
syms alpha;
syms zeta_vm zeta_hf;
for i=1:NCOUNTRIES
    s = num2str(i);
    syms(['A_',s]);
    syms(['B_',s,'_va']);
    syms(['B_',s,'_m']);
    syms(['mu_',s,'_vm']);
    for j=1:NCOUNTRIES
        s2=num2str(j);
        syms(['mu_',s,s2]);
    end
end

% final goods
syms rho_hf;
for i=1:NCOUNTRIES
    s = num2str(i);
    syms(['E_',s,'c']);
    syms(['E_',s,'i']);
    for j=1:NCOUNTRIES
        s2 = num2str(j);
        syms(['veps_',s,'c_',s2])
        syms(['veps_',s,'i_',s2])
    end
end

% other parameters
syms beta psi gama;
for i=1:NCOUNTRIES
    s = num2str(i);
    syms(['theta_',s]);
    syms(['lbar_',s]);
    syms(['delta_',s]);
end

% ---------------------------------------------------------
% equilibrium variables
% ---------------------------------------------------------

for suff_ = {'','_p'}
    suff=suff_{1};
    for i=1:NCOUNTRIES
        s = num2str(i);
        syms(['ngdp_',s,suff]);
        syms(['rgdp_',s,suff]);
        syms(['ex_',s,suff]);
        syms(['im_',s,suff]);
        syms(['nx_',s,suff]);

        syms(['y_',s,suff]);
        syms(['tfp_',s,suff]);
        syms(['va_',s,suff]);
        syms(['m_',s,suff]);
        syms(['k_',s,suff]);
        syms(['n_',s,suff]);
        syms(['k_',s,suff]);
        syms(['div_',s,suff]);
        for j=1:NCOUNTRIES
            s2=num2str(j);
            syms(['m_',s,s2,suff]);
        end
        
        syms(['fc_',s,suff]);
        syms(['fi_',s,suff]);
        for j=1:NCOUNTRIES
            s2=num2str(j);
            syms(['fc_',s,s2,suff]);
            syms(['fi_',s,s2,suff]);
        end

        syms(['w_',s,suff]);
        syms(['rk_',s,suff]);
        syms(['pfi_',s,suff]);
        syms(['pfc_',s,suff]);
        syms(['py_',s,suff]);
        syms(['q_',s,suff]);
        syms(['qm_',s,suff]);
        syms(['WEALTH_',s,suff]);
        syms(['R_',s,suff]);
        syms(['wn_',s,suff]);
        syms(['qd_',s,suff]);
    end


    for i=1:(NCOUNTRIES-1)
        s=num2str(i);
        syms(['Rx_',s,num2str(NCOUNTRIES),suff]);
        syms(['cD_',s,num2str(NCOUNTRIES),suff]);
    end

    for i=2:NCOUNTRIES
        s=num2str(i);
        syms(['rer_1',s,suff]);
        syms(['rerm_1',s,suff]);
    end
end

muc = @(c,n) c^(-gama);
mun = @(c,n,theta) -theta * n^psi;

Rs = eval(['R_',num2str(NCOUNTRIES)]);
Rs_p = eval(['R_',num2str(NCOUNTRIES),'_p']);

% ---------------------------------------------------------
% create vectors of states and controls
% ---------------------------------------------------------

% nx = 1 x (NCOUNTRIES-1) + 2 x NCOUNTRIES
% ny = 14 x NCOUNTRIES + 3 x NCOUNTRIES^2 + 3 x (NCOUNTRIES-1)
%
% TOTAL ENDOGENOUS VARIABLES:
% 16 x NCOUNTRIES + 3 x NCOUNTRIES^2 + 4 x (NCOUNTRIES-1)

% exo states, t
% nz = 3 x NCOUNTRIES)
ZZ=[];
for i=1:NCOUNTRIES
    s=num2str(i);
    ZZ = [ZZ,eval(['tfp_',s])];
end
nz=length(ZZ);

% exo states, t+1
ZZp=[];
indsZ = struct;
for i=1:nz
    name = char(ZZ(i));
    ZZp= [ZZp,eval([name,'_p'])];
    indsZ.(char(ZZ(i)))=i;
end

% endo states, t: (k_i, qm_i) for all countries, 
% WEALTH_i for all i<NCOUNTRIES, rerm_1i for all i>1
% nx = 1 x (NCOUNTRIES-1) + 2 x NCOUNTRIES
XX=[];
for i=1:NCOUNTRIES
    s=num2str(i);
    if i<NCOUNTRIES
        XX = [XX, eval(['WEALTH_',s])];
    end
    XX=[ XX, eval(['k_',s]), eval(['qm_',s])];
    %    if i>1
    %    XX = [XX, eval(['rerm_1',s])];
    %end
end
nx = length(XX);

% endo states, t+1
XXp=[];
indsX = struct;
for i=1:nx
    name = char(XX(i));
    XXp= [XXp,eval([name,'_p'])];
    indsX.(char(XX(i)))=i;
end

% controls, t
% ny = 14 x NCOUNTRIES + 3 x NCOUNTRIES^2 + 3 x (NCOUNTRIES-1)
YY=[];
for i=1:NCOUNTRIES
    s=num2str(i);
    YY=[ YY, ... 
         eval(['y_',s]), ...
         eval(['va_',s]), ...
         eval(['n_',s]) ,...
         eval(['m_',s]) ,...
         eval(['div_',s]) ,...
         eval(['fc_',s]) ,...
         eval(['fi_',s]) ,...
         eval(['pfi_',s]) ,...
         eval(['pfc_',s]) ,...
         eval(['py_',s]) ,...
         eval(['rk_',s]) ,...
         eval(['w_',s]) ,...
         eval(['R_',s]) ,...
         eval(['q_',s]),...
         eval(['wn_',s]),...
         eval(['qd_',s])];

    for j=1:NCOUNTRIES
        s2=num2str(j);
        YY=[ YY, ...
             eval(['m_',s,s2]) ...
             eval(['fc_',s,s2]) ...
             eval(['fi_',s,s2])];
    end

    if i<NCOUNTRIES
        YY=[ YY, ...
             eval(['Rx_',s,num2str(NCOUNTRIES)]) ...
             eval(['cD_',s,num2str(NCOUNTRIES)]) ];
    end
    if i>=2
        YY=[YY, eval(['rer_1',s])];
    end
end
ny = length(YY);

% controls, t+1
YYp=[];
indsY = struct;
for i=1:ny
    name = char(YY(i));
    YYp= [YYp,eval([name,'_p'])];
    indsY.(char(YY(i)))=i;
end

% *****************************************************************************
% declare analytical equilibrium conditions
% *****************************************************************************

FF=[];
nf=0;
indsF=struct;

% gross output producion functions: 3 x NCOUNTRIES
% retail output: 2 x NCOUNTRIES
% intratemporal marginal product conditions for producers: NCOUNTRIES + 3*NCOUNTRIES^2
% some other stuff: 6 x NCOUNTRIES
% real exchange rates: NCOUNTRIES-1
% budget constraints: NCOUNTRIES
% stock pricing and stock returns: 2 x NCOUNTRIES
% combined euler equations: NCOUNTRIES-1
% placeholders for price lags: NCOUNTRIES
% placeholders for excess returns: NCOUNTRIES-1
% placeholders for consumption differentials: NCOUNTRIES-1
% wage income and stock income: NCOUNTRIES*2
%
% TOTAL EQUATIONS:
% 16 x NCOUNTRIES + 3 x NCOUNTRIES^2 + 4 x (NCOUNTRIES-1)

% gross output production functions: 3 x NCOUNTRIES
for i=1:NCOUNTRIES
    s=num2str(i);

    % value added
    nf=nf+1;
    f = exp( eval(['va_',s]) ) - ...
        exp( eval(['tfp_',s]) ) * eval(['B_',s,'_va']) * ...
        (exp( eval(['k_',s]) ))^alpha * ...
        (exp( eval(['n_',s]) ))^(1-alpha);
    FF=[FF,f];
    indsF.(['va_',s]) = nf;

    % intermediate bundle
    nf=nf+1;
    tmp=0;
    for j=1:NCOUNTRIES
        s2 = num2str(j);
        tmp = tmp + (eval(['mu_',s,s2]))^(1-zeta_hf) * ...
              (exp( eval(['m_',s,s2]) ))^zeta_hf;
    end
    f = exp( eval(['m_',s]) ) - ...
        eval(['B_',s,'_m']) * tmp^(1/zeta_hf);
    FF=[FF,f];
    indsF.(['m_',s]) = nf;

    % gross output
    nf=nf+1;
    f = exp( eval(['y_',s]) ) - ...
        eval(['A_',s]) * (...
            (eval(['mu_',s,'_vm']))^(1-zeta_vm) * ...
            (exp( eval(['va_',s]) ))^zeta_vm + ...
            (1-eval(['mu_',s,'_vm']))^(1-zeta_vm) * ...
            (exp( eval(['m_',s]) ))^zeta_vm...
            )^(1/zeta_vm);
    FF=[FF,f];
    indsF.(['y_',s]) = nf;

end

% retail output: 2 x NCOUNTRIES
for i=1:NCOUNTRIES
    s = num2str(i);
    tmp1 = 0;
    tmp2 = 0;
    for j=1:NCOUNTRIES
        s2 = num2str(j);
        tmp1 = tmp1 + (eval(['veps_',s,'c_',s2]))^(1-rho_hf) * ...
               (exp( eval(['fc_',s,s2]) ))^rho_hf;
        tmp2 = tmp2 + (eval(['veps_',s,'i_',s2]))^(1-rho_hf) * ...
               (exp( eval(['fi_',s,s2]) ))^rho_hf;
    end

    nf=nf+1;
    f1 = exp( eval(['fc_',s]) ) - ...
        eval(['E_',s,'c']) * tmp1^(1/rho_hf);
    indsF.(['fc_',s]) = nf;

    nf=nf+1;
    f2 = exp( eval(['fi_',s]) ) - ...
        eval(['E_',s,'i']) * tmp2^(1/rho_hf);
    indsF.(['fi_',s]) = nf;

    FF = [FF,f1,f2];
end

% intratemporal marginal product conditions for producers:
% NCOUNTRIES + 3*NCOUNTRIES*NCOUNTRIES
for i=1:NCOUNTRIES
    s = num2str(i);
    
    nf=nf+1;
    f = exp( eval(['w_',s]) ) - ...
        exp( eval(['py_',s]) ) * (eval(['mu_',s,'_vm']))^(1-zeta_vm) * ...
        (1-alpha) * (eval(['A_',s]))^zeta_vm * eval(['B_',s,'_va']) * ...
       ( exp( eval(['y_',s]) ) / exp( eval(['va_',s]) ) )^(1-zeta_vm) * ...
     exp( eval(['tfp_',s]) ) * ...
    ( exp( eval(['k_',s]) ) / exp( eval(['n_',s]) ) )^alpha;
    indsF.(['mpn_',s,s2]) = nf;
    FF = [FF,f];

    for j=1:NCOUNTRIES
        nf = nf+1;
        s2 = num2str(j);
        f = exp( eval(['py_',s2]) ) - ...
            exp( eval(['py_',s]) ) * (1-eval(['mu_',s,'_vm']))^(1-zeta_vm) * ...
            eval(['mu_',s,s2])^(1-zeta_hf) * (eval(['A_',s]))^zeta_vm * ...
            (eval(['B_',s,'_m']))^zeta_hf * ...
            ( exp( eval(['y_',s]) ) / exp( eval(['m_',s]) ) )^(1-zeta_vm) * ...
            ( exp( eval(['m_',s]) ) / exp( eval(['m_',s,s2]) ))^(1-zeta_hf);
        indsF.(['mpm_',s,s2]) = nf;
        FF = [FF,f];
    end
end

for i=1:NCOUNTRIES
    s=num2str(i);
    for j=1:NCOUNTRIES
        s2 = num2str(j);
        nf = nf+1;
        f = exp( eval(['py_',s2]) ) - exp(eval(['pfc_',s])) * ...
            (eval(['veps_',s,'c_',s2]))^(1-rho_hf) * ...
            (eval(['E_',s,'c']))^rho_hf * ...
            ( exp( eval(['fc_',s]) ) / exp( eval(['fc_',s,s2]) ) )^(1-rho_hf);
        indsF.(['mcm_',s,s2]) = nf;
        FF = [FF,f];

        nf = nf+1;
        f = exp( eval(['py_',s2]) ) - exp(eval(['pfi_',s])) * ...
            (eval(['veps_',s,'i_',s2]))^(1-rho_hf) * ...
            (eval(['E_',s,'i']))^rho_hf * ...
            ( exp( eval(['fi_',s]) ) / exp( eval(['fi_',s,s2]) ) )^(1-rho_hf);
        indsF.(['mim_',s,s2]) = nf;
        FF = [FF,f];

    end
end

% some other stuff: 6 x NCOUNTRIES
for i=1:NCOUNTRIES

    s=num2str(i);

    % consumption-leisure first order conditions
    nf=nf+1;
    f = exp( eval(['w_',s]) ) / exp( eval(['pfc_',s]) ) + ...
        mun( exp(eval(['fc_',s])),(exp(eval(['n_',s]))/eval(['lbar_',s])),eval(['theta_',s]) ) / ...
        muc( exp(eval(['fc_',s])),(exp(eval(['n_',s]))/eval(['lbar_',s])) );
    indsF.(['mucmun_',s])=nf;
    FF=[FF,f];

    % laws of motion for capital
    nf = nf+1;
    f = exp( eval(['k_',s,'_p']) ) - ...
        ((1 - eval(['delta_',s])) * exp( eval(['k_',s]) ) + ...
         exp( eval(['fi_',s]) ) );
    indsF.(['lomk_',s])=nf;
    FF=[FF,f];

    % return on capital
    nf = nf+1;
    f = exp( eval(['rk_',s]) ) - exp( eval(['py_',s]) ) * ...
        (eval(['mu_',s,'_vm']))^(1-zeta_vm) * alpha * ...
        (eval(['A_',s]))^zeta_vm * eval(['B_',s,'_va']) * ...
        ( exp( eval(['y_',s]) ) / exp( eval(['va_',s]) ) )^(1-zeta_vm) ...
        * exp( eval(['tfp_',s]) ) * ...
        ( exp( eval(['k_',s]) ) / exp( eval(['n_',s]) ) )^(alpha-1);
    indsF.(['rk_',s]) = nf;
    FF = [FF,f];

    % Euler equation for producer
    nf = nf+1;
    f = beta*muc(exp(eval(['fc_',s,'_p'])),exp( eval(['n_',s,'_p']))) /...
        exp(eval(['pfc_',s,'_p'])) * ...
        (exp(eval(['rk_',s,'_p'])) + exp(eval(['pfi_',s,'_p'])) * (1-eval(['delta_',s]))) - ...
        exp(eval(['pfi_',s])) * muc(exp(eval(['fc_',s])),exp(eval(['n_',s])))/...
        exp(eval(['pfc_',s]));
    %end
    indsF.(['eulerk_',s]) = nf;
    FF = [FF,f];  

    % mkt clearing
    tmp1=0;
    tmp2=0;
    for j=1:NCOUNTRIES
        s2 = num2str(j);
        tmp1 = tmp1 + exp( eval(['m_',s2,s]) ) + ...
              exp( eval(['fc_',s2,s]) ) + exp( eval(['fi_',s2,s]));   
        tmp2 = tmp2 + exp( eval(['py_',s2]) ) * exp( eval(['m_',s,s2]) );
    end
    nf=nf+1;
    f = exp( eval(['y_',s]) ) - tmp1;
    indsF.(['mktclear_',s]) = nf;
    FF=[FF,f];

    % dividends
    nf=nf+1;
    f = eval(['div_',s]) - ...
        (exp(eval(['py_',s]))*exp(eval(['y_',s])) - ...
         exp(eval(['w_',s]))*exp(eval(['n_',s])) - tmp2 - ...
         exp(eval(['pfi_',s])) * exp(eval(['fi_',s])));
    indsF.(['div_',s]) = nf;
    FF=[FF,f];
end

% real exchange rates: NCOUNTRIES-1
for i=2:NCOUNTRIES
    s=num2str(i);
    nf=nf+1;
    f = exp(eval(['rer_1',s])) - ...
        exp(eval(['pfc_',s]))/exp(pfc_1);
    indsF.(['rer_1',s]) = nf;
    FF=[FF,f];
end

% budget constraints: NCOUNTRIES
snc = num2str(NCOUNTRIES);
for i=1:(NCOUNTRIES-1)
    s = num2str(i);
    nf = nf+1;
    
    % lump-sum tax/transfer related to wedge
    f = eval(['WEALTH_',s]) * exp(eval(['R_',snc]))  + ...
        eval(['div_',s]) + exp(eval(['w_',s])) * exp(eval(['n_',s])) - ...
        exp(eval(['pfc_',s])) * exp(eval(['fc_',s])) - eval(['WEALTH_',s,'_p']);
    indsF.(['bc_',s]) = nf;
    FF = [FF,f];
end

% stock pricing and stock returns: 2 x NCOUNTRIES
for i=1:NCOUNTRIES
    s = num2str(i);

    f1 = beta * (muc(exp(eval(['fc_',snc,'_p'])),exp(eval(['n_',snc,'_p']))) / exp(eval(['pfc_',snc,'_p']))) * ...
         (exp(eval(['q_',s,'_p'])) + eval(['div_',s,'_p'])) - ...
         exp(eval(['q_',s])) * ...
         (muc(exp(eval(['fc_',snc])),exp(eval(['n_',snc]))) / exp(eval(['pfc_',snc])));

    f2 = exp(eval(['q_',s])) + eval(['div_',s]) - ...
         exp(eval(['qm_',s]))*exp(eval(['R_',s]));

    nf = nf+1;
    FF=[FF,f1];
    indsF.(['q_',s]) = nf;

    nf = nf+1;
    FF=[FF,f2];
    indsF.(['R_',s]) = nf;
end

% combined euler equations: NCOUNTRIES-1
tmp = beta*exp(eval(['R_',snc','_p'])) * ...
      muc(exp(eval(['fc_',snc,'_p'])),exp(eval(['n_',snc,'_p']))) / ...
      exp(eval(['pfc_',snc,'_p'])) - ...
      muc(exp(eval(['fc_',snc])),exp(eval(['n_',snc]))) / ...
      exp(eval(['pfc_',snc]));

for i=1:(NCOUNTRIES-1)
    s = num2str(i);
    nf = nf+1;
    if(i==1)
        f = (beta*exp( eval(['R_',snc,'_p']) )*muc(exp(fc_1_p),exp(n_1_p)) - ...
             muc(exp(fc_1),exp(n_1))) - tmp;
    else
        f = beta*exp(eval(['R_',snc,'_p'])) * ...
            muc(exp(eval(['fc_',s,'_p'])),exp(eval(['n_',s,'_p']))) / ...
            exp(eval(['pfc_',s,'_p'])) - ...
            muc(exp(eval(['fc_',s])),exp(eval(['n_',s]))) / ...
            exp(eval(['pfc_',s])) - tmp;
    end
    indsF.(['eecomb_',s]) = nf;
    FF = [FF,f];
end

% numeraire
nf=nf+1;
f = exp(pfc_1) - 1.0;
FF=[FF,f];
indsF.('numeraire') = nf;

% placeholders for price lags: NCOUNTRIES
for i=1:NCOUNTRIES
    s = num2str(i);
    nf = nf+1;
    f = 0;
    FF = [FF,f];
    indsF.(['qm_',s]) = nf;
end

% placeholders for excess returns: NCOUNTRIES-1
for i=1:(NCOUNTRIES-1)
    s = num2str(i);
    nf = nf+1;
    %f = eval(['Rx_',s,snc]) - (exp(eval(['R_',s])) - exp(eval(['R_',snc])));
    f = 0;
    FF = [FF,f];
    indsF.(['Rx_',s,snc]) = nf;
end

% placeholders for consumption differentials: NCOUNTRIES-1
for i=1:(NCOUNTRIES-1)
    s = num2str(i);
    nf = nf+1;
    f = 0;
    FF = [FF,f];
    indsF.(['cD_',s,snc]) = nf;
end

% wage income and stock income
for i=1:(NCOUNTRIES)
    s=num2str(i);

    nf=nf+1;
    f=exp(eval(['wn_',s])) - exp(eval(['w_',s])) * exp(eval(['n_',s]));
    FF=[FF,f];
    indsF.(['wn_',s]) = nf;

    nf=nf+1;
    f=eval(['qd_',s]) - (exp(eval(['q_',s])) + eval(['div_',s]));
    FF=[FF,f];
    indsF.(['qd_',s]) = nf;
end
