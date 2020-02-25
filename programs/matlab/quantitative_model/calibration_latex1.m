fid = fopen([output_path,'calibration1.tex'],'wb');

fprintf(fid,'\\begin{table}[h!]\n');
fprintf(fid,'\\renewcommand{\\arraystretch}{1.2}\n');
fprintf(fid,'\\begin{center}\n');
fprintf(fid,['\\caption{Parameters fixed across calibrations}\n']);
fprintf(fid,['\\label{tab:cal1}\n']);
%fprintf(fid,'\\small\n');
fprintf(fid,'\\begin{tabular}{lc}\n');
fprintf(fid,'\\toprule\n');
fprintf(fid,'Parameter & Value\\\\\n');
fprintf(fid,'\\midrule\n');
fprintf(fid,'$\\beta$ & %0.2f\\\\\n',beta);
fprintf(fid,'$\\eta$ & %0.2f\\\\\n',1/(1-zeta_vm));
fprintf(fid,'$\\zeta$ & %0.2f\\\\\n',1.00);
fprintf(fid,'$\\rho$ & %0.2f\\\\\n',1.00);
fprintf(fid,'$\\delta$ & %0.2f\\\\\n',delta_1);
fprintf(fid,'$\\alpha$ & %0.2f\\\\\n',alpha);
fprintf(fid,'$\\gamma$ & %0.2f\\\\\n',1.00);
fprintf(fid,'$\\varphi$ & %0.2f\\\\\n',1.00);
fprintf(fid,'$');
fprintf(fid,'\\tau$ & $\\begin{bmatrix}%0.7f &%0.7f &%0.7f &%0.7f\\end{bmatrix}$\\\\[2ex]\n',TAU1(1),TAU1(2),TAU1(3),TAU1(4));
fprintf(fid,'$P$ & $\\begin{bmatrix}');
for i=1:NCOUNTRIES
    for j=1:NCOUNTRIES
        fprintf(fid,'%0.2f',N(i,j));
        if j<NCOUNTRIES
            fprintf(fid,'&');
        elseif j==NCOUNTRIES && i<NCOUNTRIES
            fprintf(fid,'\\\\');
        end
    end
end
fprintf(fid,'\\end{bmatrix}$\\\\[6ex]\n');

fprintf(fid,'$\\Sigma$ & $\\begin{bmatrix}');
for i=1:NCOUNTRIES
    for j=1:NCOUNTRIES
        fprintf(fid,'%0.5f',SIGMA(i,j));
        if j<NCOUNTRIES
            fprintf(fid,'&');
        elseif j==NCOUNTRIES && i<NCOUNTRIES
            fprintf(fid,'\\\\');
        end
    end
end
fprintf(fid,'\\end{bmatrix}$\\\\[6ex]\n');
fprintf(fid,'\\bottomrule\n');
fprintf(fid,'\\end{tabular}\n');
%fprintf(fid,'\\normalsize\n');
fprintf(fid,'\\end{center}\n');
fprintf(fid,'\\end{table}\n');
fclose(fid);









