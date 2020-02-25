% *****************************************************************************
% assigned parameters
% *****************************************************************************

% set some nontargeted parameters

delta_fixed = 0.06;
alpha=0.36; % capital share in value added
r_base = 0.04; % interest rate
beta = 1/(1+r_base); % discount factor
gama = 1; % relative risk aversion/EIS
frisch = 1; % frisch elascicity of labor supply
psi = 1/frisch; % labor supply exponent

% elasticities of substitution
el_vm=0.05; % value added vs. intermediates
el_I_hf = 1.01; % home vs. foreign in intermediates
el_F_hf = 1.01; % home vs. foreign in final goods

if sens==1
   el_I_hf = 1.25;
   el_F_hf = 1.25;
elseif sens==2
    el_F_hf=0.75;
    el_F_hf=0.75;
elseif sens==3
    el_F_hf = 1.25;
    el_I_hf = (1.0 - 0.35*el_F_hf) / 0.65;
elseif sens==4
    el_F_hf=0.75;
    el_I_hf = (1.0 - 0.35*el_F_hf) / 0.65;
elseif sens==5 || sens==7
    el_vm = 0.99;
elseif sens==8
    gama=2;
end;

% set elasticity parameters
zeta_vm = 1-1/el_vm;
zeta_hf = 1-1/el_I_hf;
rho_hf = 1-1/el_F_hf;

% *****************************************************************************
% read text files with WIOD data and process them,
% declare empty parameter arrays
% *****************************************************************************
 
% scale everything so that US GDP = 100 (assumes that US is region 1)
usgdp = va_base(1,3);
inin_base = inin_base/usgdp;
y_base = y_base/usgdp;
va_base = va_base/usgdp;
ff_base = ff_base/usgdp;

% create arrays to store numerical parameter values and intermediate
% computations
M_base = zeros(NCOUNTRIES,1);
Bv = zeros(NCOUNTRIES,1);
Bm = zeros(NCOUNTRIES,1);
mu_ = zeros(NCOUNTRIES,NCOUNTRIES);

f_base = zeros(NCOUNTRIES,2);
D = zeros(NCOUNTRIES,2);
eps_ = zeros(NCOUNTRIES,2,NCOUNTRIES);

C_base = zeros(1,NCOUNTRIES);
X_base = zeros(1,NCOUNTRIES);
L_base = zeros(1,NCOUNTRIES);
K_base = zeros(1,NCOUNTRIES);
GDP_base = zeros(1,NCOUNTRIES);
EX_base = zeros(1,NCOUNTRIES);
IM_base = zeros(1,NCOUNTRIES);
IM_F_base = zeros(1,NCOUNTRIES);
IM_I_base = zeros(1,NCOUNTRIES);
NX_base = zeros(1,NCOUNTRIES);

Lbar = zeros(1,NCOUNTRIES);
theta = zeros(1,NCOUNTRIES);

% *****************************************************************************
% calibrate gross output production parameters, taking elasticities as given
% *****************************************************************************

% set total returns to capital as fraction of GDP
delta_ = zeros(NCOUNTRIES,1);
for i=1:NCOUNTRIES
    rdky = alpha * va_base(i,3);
        
    % set depreciation rate exogenously
    delta_(i) = delta_fixed;

    % calculate capital stock
    K_base_init = rdky/(r_base+delta_(i));

    % reshuffle consumption and investment so as to make
    % investment equal to delta*K in steady state
    dky = delta_(i)*K_base_init;
    tmpc = ff_base(i,1,:);
    tmpx = ff_base(i,2,:);
    frac = dky/sum(sum(tmpx));
    tmp = tmpx * (1-frac);
    ff_base(i,2,:) = frac*tmpx;
    ff_base(i,1,:) = ff_base(i,1,:) + (1-frac)*tmpx;
    if( abs(sum(sum(ff_base(i,1,:)+ff_base(i,2,:)))- ...
            sum(sum(tmpx+tmpc))) > TINY)
        error('Error reshuffling final demand!');
    end

    % now set sector level capital stocks
    va_base(i,1) = va_base(i,3) .* alpha ./ (r_base+delta_(i));

    % set labor
    va_base(i,2) = va_base(i,3) .* (1-alpha);
    Bv(i) = va_base(i,3) ./ ( va_base(i,1).^alpha .* va_base(i,2).^(1-alpha) );
end

% ces shares and scaling factors
for i=1:NCOUNTRIES       
    % logical array of length NCOUNTRIES, set to TRUE where
    % inputs from that (sector,country) pair are non-zero,
    % and FALSE otherwise
    mask1 = inin_base(i,:) > TINY;

    % grab index of first source country with positive inputs
    [junk,idx] = max(inin_base(i,find(mask1)));

    % calculate vector of constants (see equation 17 in calinbration notes)
    % using idx as index of reference country
    tmp1 = (inin_base(i,:)./inin_base(i,idx));
        
    % solve equation 18 for share parameter of reference country 
    mu_(i,idx) = 1/sum(tmp1(find(mask1)));

    % now set shares for other countries with positive inputs
    cnt=0;
    for j=1:NCOUNTRIES
        if(mask1(j))
            cnt = cnt+1;
            if(j ~= idx)
                if(cnt<sum(mask1))
                    mu_(i,j) = mu_(i,idx) * tmp1(j);
                else
                    mu_(i,j) = 1-sum(mu_(i,find(mask1)));
                end
            end
        end
    end

    % calculate bundle of inputs from sector r
    tmp1 = inin_base(i,:).^zeta_hf;
    tmp1(find(mask1 == false)) = 0;
    M_base(i) = (dot(mu_(i,:).^(1-zeta_hf),tmp1))^(1/zeta_hf);

    % calculate scaling factor
    tmp = sum(inin_base(i,:));
    Bm(i) = tmp/M_base(i);
    M_base(i) = Bm(i) * M_base(i);
end

% value added shares
for i=1:NCOUNTRIES
    if(noio==1)
        mu_vm(i)=1;
    else
        tmp = ((1-alpha).*Bv(i))./ ...
              (mu_(i,j))^(1-zeta_hf) ...
              .*(Bm(i).^zeta_hf)...
              .* (M_base(i)./va_base(i,3)).^(1-zeta_vm) ...
              .* (va_base(i,1)./va_base(i,2)).^(alpha) ...
              .* (M_base(i)./inin_base(i,j)).^(zeta_hf-1);
        tmp = tmp^(1/(1-zeta_vm));
        mu_vm(i) = 1/(1+tmp);
    end
end

% scaling factors
for i=1:NCOUNTRIES
    if noio==1
        A(i)=1;
    else
        A(i) = y_base(i) ./ ...
                 (mu_vm(i).^(1-zeta_vm).*va_base(i,3).^zeta_vm ...
                  +(1-mu_vm(i)).^(1-zeta_vm).*M_base(i).^zeta_vm).^(1/zeta_vm);
    end
end

% check first order conditions
for i=1:NCOUNTRIES
    
    % reconstruct base-period aggregates using calibrated parameters
    M = 0;
    for j=1:NCOUNTRIES
        %if(mu_(i,j)>1.0e-15)
        %M = M + mu_(i,j)^(1-zeta_hf) * inin_base(i,j)^zeta_hf;
            %end
        M = Bm(i) * (dot(mu_(i,:).^(1-zeta_hf),...
                                 inin_base(i,:).^(zeta_hf)))^(1/(zeta_hf));
        va = Bv(i) * (va_base(i,1)^alpha) * (va_base(i,2)^(1-alpha));
        
        if noio==1
            y=va;
        else
            y =  A(i) * (mu_vm(i)^(1-zeta_vm)*va^zeta_vm + ...
                           (1-mu_vm(i))^(1-zeta_vm)*M^zeta_vm)^(1/zeta_vm);
        end

        % zero profit
        %tmp = y - (r_base+delta_(i))*va_base(i,1) - va_base(i,2) - ...
        %      sum(sum(inin_base(i,:)));

        %if(abs(tmp)>1.0e-2)
        %    error(['Zero profit condition 1 violated at i = ',...
        %           num2str(i),'! Equilibrium condition = ',num2str(tmp)]);
        %end
        
        % FOC for labor
        tmp = mu_vm(i)^(1-zeta_vm) * (1-alpha)* A(i)^zeta_vm * Bv(i) ...
              * (y / va)^(1-zeta_vm) ...
              * (va_base(i,1)/va_base(i,2))^(alpha) - 1;
        
        if(abs(tmp)>1.0e-6)
            error(['FOC for labor 1 violated at i = ',...
                   num2str(i),'! Equilibrium condition = ',num2str(tmp)]);
        end

        % FOC for capital
        tmp = mu_vm(i)^(1-zeta_vm) * alpha * A(i)^zeta_vm * Bv(i) ...
              *(y/va)^(1-zeta_vm)*(va_base(i,1)/va_base(i,2))^(alpha-1) ...
              - (r_base+delta_(i));
        
        if(abs(tmp)>1.0e-6)
            error(['FOC for capital 1 violated at i = ',...
                   num2str(i),'! Equilibrium condition = ',num2str(tmp)]);
        end
        
        % FOCs for intermediates
        for j=1:NCOUNTRIES
            tmp = (1-mu_vm(i))^(1-zeta_vm) * ...
                  mu_(i,j)^(1-zeta_hf) ...
                  * A(i)^zeta_vm * Bm(i)^zeta_hf ...
                  * (y / M)^(1-zeta_vm) ...
                  * (M / inin_base(i,j))^(1-zeta_hf) - 1;

            if(abs(tmp)>3.0e-6)
                error(['FOC for intermediates 1 violated at (i,j) = (',...
                       num2str(i),',',num2str(j),...
                       ')! Equilibrium condition = ',num2str(tmp),'\n']);
            end
        end
    end
end
   
% *****************************************************************************
% calibrate retail parameters, taking elasticities as given
% *****************************************************************************

for i=1:NCOUNTRIES
    for f=1:2
            
        % logical array of length NCOUNTRIES, set to TRUE where
        % final use of that (sector,country) pair is non-zero,
        % and FALSE otherwise
        mask1 = ff_base(i,f,:) > TINY;

        % grab index of first source country with positive inputs
        [junk,idx] = max(ff_base(i,f,find(mask1)));

        % calculate vector of J-constants using idx as index of reference country
        tmp1 = (ff_base(i,f,:)./ff_base(i,f,idx));
        
        % solve for share parameter of reference country 
        eps_(i,f,idx) = 1/sum(tmp1(find(mask1)));

        % now set shares for other countries with positive inputs
        cnt=0;
        for j=1:NCOUNTRIES
            if(mask1(j))
                cnt = cnt+1;
                if(j ~= idx)
                    if(cnt<sum(mask1))
                        eps_(i,f,j) = eps_(i,f,idx) * tmp1(j);
                    else
                        eps_(i,f,j) = 1-sum(eps_(i,f,find(mask1)));
                    end
                end
            end
        end

        % calculate bundle of inputs from sector s
        tmp1 = ff_base(i,f,:).^rho_hf;
        tmp1(find(mask1 == false)) = 0;
        f_base(i,f) = (dot(eps_(i,f,:).^(1-rho_hf),tmp1))^(1/rho_hf);

        % calculate scaling factor
        tmp = sum(ff_base(i,f,:));
        D(i,f) = tmp/f_base(i,f);
        f_base(i,f) = D(i,f) * f_base(i,f);
    end
end

% check that quantities satisfy equilibrium conditions given calibrated parameters
for i=1:NCOUNTRIES
    for f=1:2
        
        % reconstruct base-period aggregates using calibrated parameters
        yf=0;
        for j=1:NCOUNTRIES
            if(eps_(i,f,j)>1.0e-15)
                yf = yf + eps_(i,f,j)^(1-rho_hf) * (ff_base(i,f,j)^rho_hf);
            end
        end
        yf = D(i,f) * (yf^(1/rho_hf));

        % FOCS for final use
        for j=1:NCOUNTRIES
            tmp = eps_(i,f,j)^(1-rho_hf) ...
                  * D(i,f)^(rho_hf) ...
                  * (yf / ff_base(i,f,j))^(1-rho_hf) - 1;

            if(abs(tmp)>5.0e-5)
                error(['FOC for final use 1 violated at (i,f,j) = (',...
                       num2str(i),',',num2str(f),',',num2str(j),','...
                       ')! Equilibrium condition = ',num2str(tmp)]);
            end
        end
    end
end

% *****************************************************************************
% market clearing check
% *****************************************************************************

% check national accounts in raw data
for i=1:NCOUNTRIES
    C_base(i) = f_base(i,1);
    X_base(i) = f_base(i,2);
    L_base(i) = va_base(i,2);
    K_base(i) = va_base(i,1);
    GDP_base(i) = va_base(i,3);

    EX_base(i) = 0;
    IM_base(i) = 0;
    IM_I_base(i) = 0.0;
    IM_F_base(i) = 0.0;
    for j=1:NCOUNTRIES
        if(j~=i)
            EX_base(i) = EX_base(i) + sum(ff_base(j,:,i)) + inin_base(j,i);
            IM_base(i) = IM_base(i) + sum(ff_base(i,:,j)) + inin_base(i,j);
            IM_I_base(i) = IM_I_base(i) + inin_base(i,j);
            IM_F_base(i) = IM_F_base(i) + sum(ff_base(i,:,j));
        end
    end
    NX_base(i) = EX_base(i)-IM_base(i);

    tmp = (GDP_base(i) - (C_base(i) + X_base(i) + NX_base(i))) / GDP_base(i);
    if(abs(tmp) > 1.0e-6)
        error(['Sum of value added !+ C+I+G+NX for i = ',num2str(i),...
               '! Difference = ',num2str(tmp)]);          
    end
end

% now check using production functions
for i=1:NCOUNTRIES

    % reconstruct base-period gross output (supply) using calibrated parameters
    M = 0;
    for j=1:NCOUNTRIES
        if(mu_(i,j)>1.0e-15)
            M = M + mu_(i,j)^(1-zeta_hf) * inin_base(i,j)^zeta_hf;
        end
    end
    M = Bm(i) * M^(1/zeta_hf);

    va = Bv(i) * (va_base(i,1)^alpha) * (va_base(i,2)^(1-alpha));
    y =  A(i) * (mu_vm(i)^(1-zeta_vm)*va^zeta_vm + ...
                  (1-mu_vm(i))^(1-zeta_vm)*M^zeta_vm)^(1/zeta_vm);

    % calculate base period demand
    md=0;
    fd=0;
    for j=1:NCOUNTRIES
        md = md + inin_base(j,i);
        for f=1:2
            fd = fd + ff_base(j,f,i);
        end
    end
    yd=md+fd; 

    tmp = (yd-y)/y;
    if(abs(tmp)>1.0e-7)
        error(['Market clearing violated at i = ',num2str(i)...
               '! Excess demand = ',num2str(tmp)]);
    end


end

% *****************************************************************************
% calibrate household parameters and check intratemporal FOC
% *****************************************************************************

MUc = @(C) ( (C>0.0001).*(C.^(-gama)) + (C<=0.0001).*((0.1./C).*1.0e4.^(-gama)));
Lbar = 3*L_base;
theta=zeros(1,NCOUNTRIES);

for i=1:NCOUNTRIES

    % theta(i) = MUc(C_base(i))/(L_base(i)^psi);
    theta(i) = MUc(C_base(i))/((1.0/3.0)^psi);

    % check FOC
    %tmp = MUc(C_base(i))/(theta(i)*L_base(i)^psi) - 1.0;
    tmp = MUc(C_base(i))/(theta(i)*(1.0/3.0)^psi) - 1.0;
    if(abs(tmp)>1.0e-8)
        error(['Intratemporal FOC violated at i = ',num2str(i),...
              '! Equilibrium condition = ',num2str(tmp)]);
    end
end

% *****************************************************************************
% store calibrated parameters in symbolic variables
% *****************************************************************************

for i=1:NCOUNTRIES
    s = num2str(i);
    assignin('caller',['mu_',s,'_','vm'],mu_vm(i));
    for j=1:NCOUNTRIES
        s2 = num2str(j);
        assignin('caller',['mu_',s,s2],mu_(i,j));
        assignin('caller',['veps_',s,'c','_',s2],eps_(i,1,j));
        assignin('caller',['veps_',s,'i','_',s2],eps_(i,2,j));
    end
    assignin('caller',['B_',s,'_','m'],Bm(i));
    assignin('caller',['B_',s,'_','va'],Bv(i));
    assignin('caller',['A_',s],A(i));
    assignin('caller',['E_',s,'c'],D(i,1));
    assignin('caller',['E_',s,'i'],D(i,2));
    assignin('caller',['theta_',s],theta(i));
    assignin('caller',['delta_',s],delta_(i));
    assignin('caller',['lbar_',s],Lbar(i));
end

% *****************************************************************************
% calibrate initial bondholdings
% *****************************************************************************

% assume that the input-output table data represents a steady state
% then current account = zero, which implies nx = -rb

b_base = zeros(NCOUNTRIES,1);
for i=1:NCOUNTRIES
    b_base(i) = -NX_base(i)/r_base;
end

% *****************************************************************************
% check dividends are positive
% *****************************************************************************
div_chk = y_base - L_base' - M_base - f_base(:,2);
if(sum(div_chk<0)>0)
    error('Negative dividends detected!');
end

% *****************************************************************************
% write share parameters to file
% *****************************************************************************

if sens==0
    fid = fopen([output_path,'shareparams',str,mname,'.txt'],'wb');
    fprintf(fid,'%f ',mu_);
    fprintf(fid,'\n');
    fprintf(fid,'%f ',eps_);
    fprintf(fid,'\n');
    fprintf(fid,'%f ',mu_vm);
    fprintf(fid,'\n');
    fprintf(fid,'%f ',Bm);
    fprintf(fid,'\n');
    fprintf(fid,'%f ',Bv);
    fprintf(fid,'\n');
    fprintf(fid,'%f ',A);
    fprintf(fid,'\n');
    fprintf(fid,'%f ',D);
    fprintf(fid,'\n');
    fprintf(fid,'%f ',theta);
    fprintf(fid,'\n');
    fprintf(fid,'%f ',Lbar);
    fprintf(fid,'\n');
    fclose(fid);
end

%if sens==0
%    calibration_latex;
%end
