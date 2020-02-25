%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% benchmark results (balanced-trade version)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---------------------------------------
% main table 1: changes in diversification in data, model, and sensitivity analyses

% open data files constructed in python
da = importdata('../../python/results/div_for_matlab.csv');

% open tex file
fid = fopen([output_path,'results-bal-div-base-sens.tex'],'wb');

% header info
fprintf(fid,'\\begin{table}[p]\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{threeparttable}\n');
fprintf(fid,['\\caption{Changes in diversification: data, baseline model, and sensitivity analyses}\n']);
fprintf(fid,['\\label{tab:results-bal',str,'div-sens}\n']);
%fprintf(fid,'\\small\n');

fprintf(fid,'\\begin{tabular}{l c c c c c}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'& & & \\multicolumn{3}{c}{Counterfactuals}\\\\\n');
fprintf(fid,'\\cmidrule(rl){4-6}\n');

% columns names
fprintf(fid,['Country & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering Data} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering Benchmark} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering 1. Size} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering 2. Openness} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering 3. Int. trade}\\\\\n']);
%fprintf(fid,['Country & Data & Benchmark & 1. Size & 2. Openness & 3. Int. trade\\\\\n']);

% captions for the table
captions = {'(b) Higher Armington elasticities',...
            '(c) Lower Armington elasticities',...
            '(d) Higher final Armington elasticity',...
            '(e) Higher intermediate Armington elasticity',...
            '(b) Uncorrelated shocks',...
            '(c) Cobb-Douglas production',...
            '(d) Uncorrelated shocks + Cobb-Douglas',...
            '(f) Higher risk aversion'};

%for s=[0,1,2,3,4,5,6,7,8]
for s=[0,1,2,3,4,8]

    fprintf(fid,'\\midrule\n');

    if(s==0)
        fprintf(fid,['\\multicolumn{6}{l}{\\textit{(a) Baseline}}\\\\\n']);
        a = importdata([output_path,'results',str,mnames{7},'.txt']);
        b0 = importdata([output_path,'results',str,mnames{8},'.txt']);
        b1 = importdata([output_path,'results',str,mnames{9},'.txt']);
        b2 = importdata([output_path,'results',str,mnames{10},'.txt']);
        b3 = importdata([output_path,'results',str,mnames{11},'.txt']);
    else
        fprintf(fid,['\\multicolumn{6}{l}{\\textit{',captions{s},'}}\\\\\n']);
        a = importdata([output_path,'results-s',num2str(s),str,mnames{7},'.txt']);
        b0 = importdata([output_path,'results-s',num2str(s),str,mnames{8},'.txt']);
        b1 = importdata([output_path,'results-s',num2str(s),str,mnames{9},'.txt']);
        b2 = importdata([output_path,'results-s',num2str(s),str,mnames{10},'.txt']);
        b3 = importdata([output_path,'results-s',num2str(s),str,mnames{11},'.txt']);
    end

    % data
    for i=1:NCOUNTRIES
        
        % model
        a2=100*(1-a(i,i));
        b02=100*(1-b0(i,i));
        b12=100*(1-b1(i,i));
        b22=100*(1-b2(i,i));
        b32=100*(1-b3(i,i));

        % country name
        fprintf(fid,cnames{i});

        % data
        da2=da.data(i,divcol);
        db2=da.data(i,divcol+3);
        fprintf(fid,'& %0.2f',db2-da2);
        fprintf(fid,'&');

        % benchmark
        fprintf(fid,'%0.2f',b02-a2);
        fprintf(fid,'&');

        % counter 1
        fprintf(fid,'%0.2f',b12-a2);
        fprintf(fid,'&');

        % counter 2
        fprintf(fid,'%0.2f',b22-a2);
        fprintf(fid,'&');

        % counter 3
        fprintf(fid,'%0.2f',b32-a2);
        fprintf(fid,'\\\\\n');
    end

end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\normalsize\n');

fprintf(fid,'\\begin{tablenotes}\n');
fprintf(fid,'\\small\n');
fprintf(fid,'\\item Note: all results reported above are in percentage points.\n');
fprintf(fid,'\\normalsize\n');
fprintf(fid,'\\end{tablenotes}\n');
fprintf(fid,'\\end{threeparttable}\n');

fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);

% ---------------------------------------
% main table 2: stochastic structure analyses

% open tex file
fid = fopen([output_path,'results-bal-div-base-sens2.tex'],'wb');

% header info
fprintf(fid,'\\begin{table}[p]\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{threeparttable}\n');
fprintf(fid,['\\caption{Changes in diversification: data, baseline model, and more sensitivity analyses}\n']);
fprintf(fid,['\\label{tab:results-bal',str,'div-sens2}\n']);
%fprintf(fid,'\\small\n');

fprintf(fid,'\\begin{tabular}{l c c c c c}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'& & & \\multicolumn{3}{c}{Counterfactuals}\\\\\n');
fprintf(fid,'\\cmidrule(rl){4-6}\n');

% columns names
fprintf(fid,['Country & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering Data} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering Benchmark} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering 1. Size} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering 2. Openness} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering 3. Int. trade}\\\\\n']);


for s=[0,5,6,7]

    fprintf(fid,'\\midrule\n');

    if(s==0)
        fprintf(fid,['\\multicolumn{6}{l}{\\textit{(a) Baseline}}\\\\\n']);
        a = importdata([output_path,'results',str,mnames{7},'.txt']);
        b0 = importdata([output_path,'results',str,mnames{8},'.txt']);
        b1 = importdata([output_path,'results',str,mnames{9},'.txt']);
        b2 = importdata([output_path,'results',str,mnames{10},'.txt']);
        b3 = importdata([output_path,'results',str,mnames{11},'.txt']);
    else
        fprintf(fid,['\\multicolumn{6}{l}{\\textit{',captions{s},'}}\\\\\n']);
        a = importdata([output_path,'results-s',num2str(s),str,mnames{7},'.txt']);
        b0 = importdata([output_path,'results-s',num2str(s),str,mnames{8},'.txt']);
        b1 = importdata([output_path,'results-s',num2str(s),str,mnames{9},'.txt']);
        b2 = importdata([output_path,'results-s',num2str(s),str,mnames{10},'.txt']);
        b3 = importdata([output_path,'results-s',num2str(s),str,mnames{11},'.txt']);
    end

    % data
    for i=1:NCOUNTRIES
        
        % model
        a2=100*(1-a(i,i));
        b02=100*(1-b0(i,i));
        b12=100*(1-b1(i,i));
        b22=100*(1-b2(i,i));
        b32=100*(1-b3(i,i));

        % country name
        fprintf(fid,cnames{i});

        % data
        da2=da.data(i,divcol);
        db2=da.data(i,divcol+3);
        fprintf(fid,'& %0.2f',db2-da2);
        fprintf(fid,'&');

        % benchmark
        fprintf(fid,'%0.2f',b02-a2);
        fprintf(fid,'&');

        % counter 1
        fprintf(fid,'%0.2f',b12-a2);
        fprintf(fid,'&');

        % counter 2
        fprintf(fid,'%0.2f',b22-a2);
        fprintf(fid,'&');

        % counter 3
        fprintf(fid,'%0.2f',b32-a2);
        fprintf(fid,'\\\\\n');
    end

end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\normalsize\n');

fprintf(fid,'\\begin{tablenotes}\n');
fprintf(fid,'\\small\n');
fprintf(fid,'\\item Note: all results reported above are in percentage points.\n');
fprintf(fid,'\\normalsize\n');
fprintf(fid,'\\end{tablenotes}\n');
fprintf(fid,'\\end{threeparttable}\n');

fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);

% ---------------------------------------
% equilibrium portfolio weights

% open file
fid = fopen([output_path,'results-bal',str,'weights.tex'],'wb');

% header info
fprintf(fid,'\\begin{table}[p]\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{threeparttable}\n');
fprintf(fid,['\\caption{Changes in bilateral portfolio shares: baseline model}\n']);
fprintf(fid,['\\label{tab:results-bal-weights}\n']);
%fprintf(fid,'\\small\n');

fprintf(fid,'\\begin{tabular}{c');
for i=1:NCOUNTRIES
    fprintf(fid,'c');
end
fprintf(fid,'}\n');
fprintf(fid,'\\toprule\n');

% columns: country names
fprintf(fid,'Home\\textbackslash Foreign');
for i=1:NCOUNTRIES
    fprintf(fid,'&');
    fprintf(fid,['\\multicolumn{1}{p{1.5cm}}{\\centering ',cnames{i},'}']);
end
fprintf(fid,'\\\\\n');

% write results
panels={'(a) Benchmark',...
        '(b) Counterfactual 1: changes in size only',...
        '(c) Counterfactual 2: changes in trade openness only',...
        '(d) Counterfactual 3: changes in intermediate share share only'};

results0 = importdata([output_path,'results',str,mnames{7},'.txt']);
for i=8:11
    mname=mnames{i};
    panel=panels{i-7};
    results = importdata([output_path,'results',str,mname,'.txt']);
    fprintf(fid,'\\midrule\n');
    fprintf(fid,['\\multicolumn{5}{l}{\\textit{',panel,'}}\\\\\n']);
    for i=1:NCOUNTRIES
        fprintf(fid,cnames{i});
        for j=1:NCOUNTRIES
            fprintf(fid,'&');
            fprintf(fid,'%0.2f',100*(results(i,j)-results0(i,j)));
        end
        fprintf(fid,'\\\\\n');
    end
end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\normalsize\n');

fprintf(fid,'\\begin{tablenotes}\n');
fprintf(fid,'\\small\n');
fprintf(fid,'\\item Note: all results reported above are in percentage points.\n');
fprintf(fid,'\\normalsize\n');
fprintf(fid,'\\end{tablenotes}\n');
fprintf(fid,'\\end{threeparttable}\n');

fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);

% ---------------------------------------
% unbalanced-trade version of main table

% open tex file
fid = fopen([output_path,'results-div-base-sens.tex'],'wb');

% header info
fprintf(fid,'\\begin{table}[p]\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,'\\begin{threeparttable}\n');
fprintf(fid,['\\caption{Changes in diversification: data and unbalanced-trade calibrations}\n']);
fprintf(fid,['\\label{tab:results',str,'div-sens}\n']);
%fprintf(fid,'\\small\n');

fprintf(fid,'\\begin{tabular}{l c c c c c c}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'& & & \\multicolumn{4}{c}{Counterfactuals}\\\\\n');
fprintf(fid,'\\cmidrule(rl){4-7}\n');

% columns names
fprintf(fid,['Country & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering Data} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering Benchmark} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering 1. Size} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering 2. Openness} & ',...
             '\\multicolumn{1}{p{2.1cm}}{\\centering 3. Int. trade} &',....
             '\\multicolumn{1}{p{2.1cm}}{\\centering 4. Imbalances}\\\\\n']);


for s=[0,5,6,7]

    fprintf(fid,'\\midrule\n');

    if(s==0)
        fprintf(fid,['\\multicolumn{6}{l}{\\textit{(a) Baseline}}\\\\\n']);
        a = importdata([output_path,'results',str,mnames{1},'.txt']);
        b0 = importdata([output_path,'results',str,mnames{2},'.txt']);
        b1 = importdata([output_path,'results',str,mnames{3},'.txt']);
        b2 = importdata([output_path,'results',str,mnames{4},'.txt']);
        b3 = importdata([output_path,'results',str,mnames{5},'.txt']);
        b4 = importdata([output_path,'results',str,mnames{6},'.txt']);
    else
        fprintf(fid,['\\multicolumn{6}{l}{\\textit{',captions{s},'}}\\\\\n']);
        a = importdata([output_path,'results-s',num2str(s),str,mnames{1},'.txt']);
        b0 = importdata([output_path,'results-s',num2str(s),str,mnames{2},'.txt']);
        b1 = importdata([output_path,'results-s',num2str(s),str,mnames{3},'.txt']);
        b2 = importdata([output_path,'results-s',num2str(s),str,mnames{4},'.txt']);
        b3 = importdata([output_path,'results-s',num2str(s),str,mnames{5},'.txt']);
        b4 = importdata([output_path,'results-s',num2str(s),str,mnames{6},'.txt']);
    end

    % data
    for i=1:NCOUNTRIES
        
        % model
        a2=100*(1-a(i,i));
        b02=100*(1-b0(i,i));
        b12=100*(1-b1(i,i));
        b22=100*(1-b2(i,i));
        b32=100*(1-b3(i,i));
        b42=100*(1-b4(i,i));

        % country name
        fprintf(fid,cnames{i});

        % data
        da2=da.data(i,divcol);
        db2=da.data(i,divcol+3);
        fprintf(fid,'& %0.2f',db2-da2);
        fprintf(fid,'&');

        % benchmark
        fprintf(fid,'%0.2f',b02-a2);
        fprintf(fid,'&');

        % counter 1
        fprintf(fid,'%0.2f',b12-a2);
        fprintf(fid,'&');

        % counter 2
        fprintf(fid,'%0.2f',b22-a2);
        fprintf(fid,'&');

        % counter 3
        fprintf(fid,'%0.2f',b32-a2);
        fprintf(fid,'&');

        % counter 4
        fprintf(fid,'%0.2f',b42-a2);
        fprintf(fid,'\\\\\n');

    end

end

fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
fprintf(fid,'\\normalsize\n');

fprintf(fid,'\\begin{tablenotes}\n');
fprintf(fid,'\\small\n');
fprintf(fid,'\\item Note: all results reported above are in percentage points.\n');
fprintf(fid,'\\normalsize\n');
fprintf(fid,'\\end{tablenotes}\n');
fprintf(fid,'\\end{threeparttable}\n');

fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{table}\n');

fclose(fid);
