% *****************************************************************************
% declare symbolic variables
% *****************************************************************************

% ---------------------------------------------------------
% parameters
% ---------------------------------------------------------

syms beta theta gama;

% ---------------------------------------------------------
% equilibrium variables
% ---------------------------------------------------------

for suff_ = {'','_p'}
    suff=suff_{1};
    for i=1:NCOUNTRIES
        s = num2str(i);
        syms(['y_',s,suff]);
        syms(['c_',s,suff]);
        syms(['WEALTH_',s,suff]);
        syms(['q_',s,suff]);
        syms(['qm_',s,suff]);
        syms(['R_',s,suff]);
        syms(['d_',s,suff]);
    end

    for j=1:(NCOUNTRIES-1)
        s=num2str(j);
        syms(['Rx_',s,num2str(NCOUNTRIES),suff]);
    end

    for j=1:(NCOUNTRIES-1)
        s=num2str(j);
        syms(['cD_',s,num2str(NCOUNTRIES),suff]);
    end
end

muc = @(c) (c^(-gama));
Rs = eval(['R_',num2str(NCOUNTRIES)]);
Rs_p = eval(['R_',num2str(NCOUNTRIES),'_p']);

% ---------------------------------------------------------
% create vectors of states and controls
% ---------------------------------------------------------

% exo states, t
ZZ=[];
for i=1:NCOUNTRIES
    s=num2str(i);
    ZZ = [ZZ,eval(['y_',s])];
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

% endo states, t: (WEALTH_i, qm_i) for all countries, rerm_1j
% for all j
XX=[];
for i=1:NCOUNTRIES-1
    s=num2str(i);
    XX=[ XX, eval(['WEALTH_',s]) ];
end
for i=1:NCOUNTRIES
    s=num2str(i);
    XX=[XX, eval(['qm_',s]) ];
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
YY=[];
for i=1:NCOUNTRIES
    s=num2str(i);
    YY=[ YY, ...
         eval(['c_',s]), ...
         eval(['q_',s]), ...
         eval(['d_',s]), ...
         eval(['R_',s])];
end
for i=1:(NCOUNTRIES-1)
    s=num2str(i);
    YY=[ YY, eval(['Rx_',s,num2str(NCOUNTRIES)])];
end
for i=1:(NCOUNTRIES-1)
    s=num2str(i);
    YY = [YY, eval(['cD_',s,num2str(NCOUNTRIES)]) ];
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

% budget constraints (xNCOUNTRIES)
for i=1:NCOUNTRIES-1
    s=num2str(i);
    nf = nf+1;
    f = exp(eval(['y_',s])) + exp(Rs) * eval(['WEALTH_',s]) - ...
        exp(eval(['c_',s])) - eval(['WEALTH_',s,'_p']);
    FF = [FF,f];
    indsF.(['bc_',s]) = nf;
end

% output market clearing (x1)
f = 0;
for i=1:NCOUNTRIES
    s=num2str(i);
    f = f + exp(eval(['y_',s])) - exp(eval(['c_',s]));
    %f = f + eval(['WEALTH_',s,'_p']);
end
nf = nf+1;
FF = [FF,f];
indsF.('mktcl') = nf;

% return definitions (xNCOUNTRIES)
for i=1:NCOUNTRIES
    s=num2str(i);
    nf = nf+1;
    f = (exp(eval(['q_',s])) + exp(eval(['d_',s]))) - ...
        exp(eval(['qm_',s])) * exp(eval(['R_',s]));
    FF = [FF,f];
    indsF.(['R_',s]) = nf;
end

% dividend definitions
for i=1:NCOUNTRIES
    s=num2str(i);
    nf = nf+1;
    f = theta*exp(eval(['y_',s])) - exp(eval(['d_',s]));
    FF = [FF,f];
    indsF.(['d_',s]) = nf;
end

% stock pricing equations (xNCOUNTRIES)
for j=1:NCOUNTRIES
    s = num2str(j);
    nf = nf+1;
    f = (exp(eval(['q_',s,'_p'])) + exp(eval(['d_',s,'_p']))) ...
        * beta * muc(exp(c_1_p)) - exp(eval(['q_',s])) * muc(exp(c_1));
    FF = [FF,f];
    indsF.(['q_',s]) = nf;
end

% combined euler equations (xNCOUNTRIES-1)
for j=2:NCOUNTRIES
    s=num2str(j);
    nf = nf+1;
    f = (beta*muc(exp(c_1_p))*exp(Rs_p) / muc(exp(c_1))) ...
        - (beta*muc(exp(eval(['c_',s,'_p'])))*exp(Rs_p) / muc(exp(eval(['c_',s]))));
    FF = [FF,f];
    indsF.(['ee_1',s]) = nf;
end

% placeholders for lagged prices
for i=1:NCOUNTRIES
    s=num2str(i);
    nf = nf+1;
    f=0;
    FF=[FF,0];
    indsF.(['qm_',s]) = nf;
end

% placeholders for excess returns
for i=1:(NCOUNTRIES-1)
    s=num2str(i);
    nf = nf+1;
    f=0;
    FF=[FF,0];
    indsF.(['Rx_',s,num2str(NCOUNTRIES)]) = nf;
end

% placeholders for consumption differentials
for i=1:(NCOUNTRIES-1)
    s=num2str(i);
    nf = nf+1;
    f=0;
    FF=[FF,0];
    indsF.(['cD_',s,num2str(NCOUNTRIES)]) = nf;
end
