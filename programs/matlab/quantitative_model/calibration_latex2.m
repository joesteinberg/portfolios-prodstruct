%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pars1995 = struct;
pars2011 = struct;
parsalt1 = struct;
parsalt2 = struct;
parsalt3 = struct;
parsalt4 = struct;

fid = fopen([output_path,'shareparams-bal-bench-1995.txt'],'rb');
pars1995.mu_ = fscanf(fid,'%f ',[NCOUNTRIES,NCOUNTRIES]);
tmp = fscanf(fid,'%f ',[NCOUNTRIES,2*NCOUNTRIES]);
pars1995.eps_ = reshape(tmp,[NCOUNTRIES,2,NCOUNTRIES]);
pars1995.mu_vm = fscanf(fid,'%f ',[NCOUNTRIES]);
pars1995.Bm = fscanf(fid,'%f ',[NCOUNTRIES]);
pars1995.Bv = fscanf(fid,'%f ',[NCOUNTRIES]);
pars1995.A = fscanf(fid,'%f ',[NCOUNTRIES]);
pars1995.D = fscanf(fid,'%f ',[NCOUNTRIES,2]);
pars1995.theta = fscanf(fid,'%f ',[NCOUNTRIES]);
pars1995.lbar = fscanf(fid,'%f ',[NCOUNTRIES]);
fclose(fid);

fid = fopen([output_path,'shareparams-bal-bench-2011.txt'],'rb');
pars2011.mu_ = fscanf(fid,'%f ',[NCOUNTRIES,NCOUNTRIES]);
tmp = fscanf(fid,'%f ',[NCOUNTRIES,2*NCOUNTRIES]);
pars2011.eps_ = reshape(tmp,[NCOUNTRIES,2,NCOUNTRIES]);
pars2011.mu_vm = fscanf(fid,'%f ',[NCOUNTRIES]);
pars2011.Bm = fscanf(fid,'%f ',[NCOUNTRIES]);
pars2011.Bv = fscanf(fid,'%f ',[NCOUNTRIES]);
pars2011.A = fscanf(fid,'%f ',[NCOUNTRIES]);
pars2011.D = fscanf(fid,'%f ',[NCOUNTRIES,2]);
pars2011.theta = fscanf(fid,'%f ',[NCOUNTRIES]);
pars2011.lbar = fscanf(fid,'%f ',[NCOUNTRIES]);
fclose(fid);

fid = fopen([output_path,'shareparams-bal-size-counter.txt'],'rb');
parsalt1.mu_ = fscanf(fid,'%f ',[NCOUNTRIES,NCOUNTRIES]);
tmp = fscanf(fid,'%f ',[NCOUNTRIES,2*NCOUNTRIES]);
parsalt1.eps_ = reshape(tmp,[NCOUNTRIES,2,NCOUNTRIES]);
parsalt1.mu_vm = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt1.Bm = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt1.Bv = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt1.A = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt1.D = fscanf(fid,'%f ',[NCOUNTRIES,2]);
parsalt1.theta = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt1.lbar = fscanf(fid,'%f ',[NCOUNTRIES]);
fclose(fid);

fid = fopen([output_path,'shareparams-bal-trd-counter.txt'],'rb');
parsalt2.mu_ = fscanf(fid,'%f ',[NCOUNTRIES,NCOUNTRIES]);
tmp = fscanf(fid,'%f ',[NCOUNTRIES,2*NCOUNTRIES]);
parsalt2.eps_ = reshape(tmp,[NCOUNTRIES,2,NCOUNTRIES]);
parsalt2.mu_vm = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt2.Bm = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt2.Bv = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt2.A = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt2.D = fscanf(fid,'%f ',[NCOUNTRIES,2]);
parsalt2.theta = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt2.lbar = fscanf(fid,'%f ',[NCOUNTRIES]);
fclose(fid);

fid = fopen([output_path,'shareparams-bal-io-counter.txt'],'rb');
parsalt3.mu_ = fscanf(fid,'%f ',[NCOUNTRIES,NCOUNTRIES]);
tmp = fscanf(fid,'%f ',[NCOUNTRIES,2*NCOUNTRIES]);
parsalt3.eps_ = reshape(tmp,[NCOUNTRIES,2,NCOUNTRIES]);
parsalt3.mu_vm = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt3.Bm = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt3.Bv = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt3.A = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt3.D = fscanf(fid,'%f ',[NCOUNTRIES,2]);
parsalt3.theta = fscanf(fid,'%f ',[NCOUNTRIES]);
parsalt3.lbar = fscanf(fid,'%f ',[NCOUNTRIES]);
fclose(fid);

%fid = fopen([output_path,'shareparams-nx-counter.txt'],'rb');
%parsalt4.mu_ = fscanf(fid,'%f ',[NCOUNTRIES,NCOUNTRIES]);
%tmp = fscanf(fid,'%f ',[NCOUNTRIES,2*NCOUNTRIES]);
%parsalt4.eps_ = reshape(tmp,[NCOUNTRIES,2,NCOUNTRIES]);
%parsalt4.mu_vm = fscanf(fid,'%f ',[NCOUNTRIES]);
%parsalt4.Bm = fscanf(fid,'%f ',[NCOUNTRIES]);
%parsalt4.Bv = fscanf(fid,'%f ',[NCOUNTRIES]);
%parsalt4.A = fscanf(fid,'%f ',[NCOUNTRIES]);
%parsalt4.D = fscanf(fid,'%f ',[NCOUNTRIES,2]);
%parsalt4.theta = fscanf(fid,'%f ',[NCOUNTRIES]);
%parsalt4.lbar = fscanf(fid,'%f ',[NCOUNTRIES]);
%fclose(fid);

pars = {pars1995,pars2011,parsalt1,parsalt2,parsalt3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make the latex table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fid = fopen([output_path,'calibration2.tex'],'wb');

fprintf(fid,'\\begin{table}[p]\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,['\\caption{Parameters that vary across calibrations}\n']);
fprintf(fid,['\\label{tab:cal2}\n']);
%fprintf(fid,'\\scriptsize\n');
fprintf(fid,'\\begin{tabular}{cccccc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'& & & \\multicolumn{3}{c}{Counterfactuals}\\\\\n');
fprintf(fid,'\\cmidrule(rl){4-6}\n');

% columns names
fprintf(fid,['Parameter & ',...
             '\\multicolumn{1}{p{4.0cm}}{\\centering 1995 Benchmark} & ',...
             '\\multicolumn{1}{p{4.0cm}}{\\centering 2011 Benchmark} & ',...
             '\\multicolumn{1}{p{4.0cm}}{\\centering 1. Size} & ',...
             '\\multicolumn{1}{p{4.0cm}}{\\centering 2. Openness} & ',...
             '\\multicolumn{1}{p{4.0cm}}{\\centering 3. Int. trade}\\\\\n']);

%fprintf(fid,['Parameter & 1995 Benchmark & 2011 Benchmark &  Counterfactual ' ...
%             '1 & Counterfactual 2 & Counterfactual 3 & Counterfactual 4\\\\\n']);
fprintf(fid,'\\midrule\n');

% nbar
fprintf(fid,'$\\Theta_i$');
for p = 1:5
    fprintf(fid,'& $\\begin{bmatrix}');
    for i=1:NCOUNTRIES
        fprintf(fid,'%0.2f',pars{p}.lbar(i));
        if i<NCOUNTRIES
            fprintf(fid,'&');
        end
    end
    fprintf(fid,'\\end{bmatrix}$');
end
fprintf(fid,'\\\\[3ex]\n');

% theta
fprintf(fid,'$\\theta_i$');
for p = 1:5
    fprintf(fid,'& $\\begin{bmatrix}');
    for i=1:NCOUNTRIES
        fprintf(fid,'%0.2f',pars{p}.theta(i));
        if i<NCOUNTRIES
            fprintf(fid,'&');
        end
    end
    fprintf(fid,'\\end{bmatrix}$');
end
fprintf(fid,'\\\\[3ex]\n');

% Y_i
%fprintf(fid,'$Y_i$');
%for p = 1:6
%    fprintf(fid,'& $\\begin{bmatrix}');
%    for i=1:NCOUNTRIES
%        fprintf(fid,'%0.2f',pars{p}.A(i));
%        if i<NCOUNTRIES
%            fprintf(fid,'&');
%        end
%    end
%    fprintf(fid,'\\end{bmatrix}$');
%end
%fprintf(fid,'\\\\[2ex]\n');

% V_i
%fprintf(fid,'$V_i$');
%for p = 1:6
%    fprintf(fid,'& $\\begin{bmatrix}');
%    for i=1:NCOUNTRIES
%        fprintf(fid,'%0.2f',pars{p}.Bv(i));
%        if i<NCOUNTRIES
%            fprintf(fid,'&');
%        end
%    end
%    fprintf(fid,'\\end{bmatrix}$');
%end
%fprintf(fid,'\\\\[3ex]\n');

% M_i
%fprintf(fid,'$M_i$');
%for p = 1:6
%    fprintf(fid,' & $\\begin{bmatrix}');
%    for i=1:NCOUNTRIES
%        fprintf(fid,'%0.2f',pars{p}.Bm(i));
%        if i<NCOUNTRIES
%            fprintf(fid,'&');
%        end
%    end
%    fprintf(fid,'\\end{bmatrix}$');
%end
%fprintf(fid,'\\\\[2ex]\n');

% C_i
%fprintf(fid,'$C_i$');
%for p = 1:6
%    fprintf(fid,' & $\\begin{bmatrix}');
%    for i=1:NCOUNTRIES
%        fprintf(fid,'%0.2f',pars{p}.D(i,1));
%        if i<NCOUNTRIES
%            fprintf(fid,'&');
%        end
%    end
%    fprintf(fid,'\\end{bmatrix}$');
%end
%fprintf(fid,'\\\\[2ex]\n');

% X_i
%fprintf(fid,'$X_i$');
%for p = 1:6
%    fprintf(fid,' & $\\begin{bmatrix}');
%    for i=1:NCOUNTRIES
%        fprintf(fid,'%0.2f',pars{p}.D(i,2));
%        if i<NCOUNTRIES
%            fprintf(fid,'&');
%        end
%    end
%    fprintf(fid,'\\end{bmatrix}$');
%end
%fprintf(fid,'\\\\[2ex]\n');

% upsilon_i
fprintf(fid,'$\\upsilon_i$');
for p = 1:5
    fprintf(fid,' & $\\begin{bmatrix}');
    for i=1:NCOUNTRIES
        fprintf(fid,'%0.2f',pars{p}.mu_vm(i));
        if i<NCOUNTRIES
            fprintf(fid,'&');
        end
    end
    fprintf(fid,'\\end{bmatrix}$');
end
fprintf(fid,'\\\\[3ex]\n');

% mu
fprintf(fid,'$\\mu_{i,j}$');
for p = 1:5
    fprintf(fid,' & $\\begin{bmatrix}');
    for i=1:NCOUNTRIES
        for j=1:NCOUNTRIES
            fprintf(fid,'%0.2f',pars{p}.mu_(i,j));
            if j<NCOUNTRIES
                fprintf(fid,'&');
            elseif j==NCOUNTRIES && i<NCOUNTRIES
                fprintf(fid,'\\\\');
            end
            
        end
    end
    fprintf(fid,'\\end{bmatrix}$');
end
fprintf(fid,'\\\\[8ex]\n');


% omega c
fprintf(fid,'$\\omega_{i,c,j}$');
for p = 1:5
    fprintf(fid,' & $\\begin{bmatrix}');
    for i=1:NCOUNTRIES
        for j=1:NCOUNTRIES
            fprintf(fid,'%0.2f',pars{p}.eps_(i,1,j));
            if j<NCOUNTRIES
                fprintf(fid,'&');
            elseif j==NCOUNTRIES && i<NCOUNTRIES
                fprintf(fid,'\\\\');
            end
            
        end
    end
    fprintf(fid,'\\end{bmatrix}$');
end
fprintf(fid,'\\\\[8ex]\n');

% omega x
fprintf(fid,'$\\omega_{i,x,j}$');
for p = 1:5
    fprintf(fid,' & $\\begin{bmatrix}');
    for i=1:NCOUNTRIES
        for j=1:NCOUNTRIES
            fprintf(fid,'%0.2f',pars{p}.eps_(i,2,j));
            if j<NCOUNTRIES
                fprintf(fid,'&');
            elseif j==NCOUNTRIES && i<NCOUNTRIES
                fprintf(fid,'\\\\');
            end
            
        end
    end
    fprintf(fid,'\\end{bmatrix}$');
end
fprintf(fid,'\\\\[8ex]\n');

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
%fprintf(fid,'\\normalsize\n');
fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);









