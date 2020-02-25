##################################################################################
# This script creates table 6 (GDP distribution and trade measures in benchmark
# and counterfactual input-output tables)

##################################################################################
# Imports, constants, etc.

import pandas as pd
import numpy as np
import itertools

nc=4
excel_path='../../excel/'
texpath = 'tabfig/'
countries = ['USA','ADV','EME','ROW']

names=['bench-1995',
       'bench-2011',
       'size-counter',
       'trd-counter',
       'io-counter',
       'nx-counter',
       'bal-bench-1995',
       'bal-bench-2011',
       'bal-size-counter',
       'bal-trd-counter',
       'bal-io-counter']

titles=['Benchmark 1995 IO matrix',
        'Benchmark 2011 IO matrix',
        'Counterfactual 1: change in relative size only',
        'Counterfactual 2: change in openness only',
        'Counterfactual 3: change in intermediate trade shares only',
        'Counterfactual 4: change in trade balances only',
        'Balanced-trade benchmark 1995 IO matrix',
        'Balanced-trade benchmark 2011 IO matrix',
        'Balanced-trade counterfactual 1: change in relative size only',
        'Balanced-trade counterfactual 2: change in openness only',
        'Balanced-trade counterfactual 3: change in intermediate trade shares only']

def read_iomat(fname):
    iomat = np.genfromtxt(fname=fname,dtype='float',delimiter=',',names=None)
    return iomat

def gdp_share(iomat,i):
    return 100*iomat[4,i]/sum(iomat[4,:])

def openness(iomat,i):
    trd = 0
    for j in range(0,nc):
        if j != i:
            trd = trd + iomat[i,j]+iomat[i,nc+j]+iomat[i,2*nc+j]
            trd = trd + iomat[j,i]+iomat[j,nc+i]+iomat[j,2*nc+i]
    #return 100*trd/iomat[4,i]
    return 100*trd/sum(iomat[4,:])

def itrd(iomat,i):
    trd = 0
    itrd = 0
    for j in range(0,nc):
        if j != i:
            trd = trd + iomat[i,j]+iomat[i,nc+j]+iomat[i,2*nc+j]
            trd = trd + iomat[j,i]+iomat[j,nc+i]+iomat[j,2*nc+i]
            itrd = itrd + iomat[j,i]
            itrd = itrd + iomat[i,j]
    return 100*itrd/trd

def ftrd(iomat,i):
    ff = np.sum(iomat[:,nc+i])+np.sum(iomat[:,nc*2+i])
    fm = 0
    for j in range(0,nc):
        if j != i:
            fm = fm + iomat[j,nc+i] + iomat[j,nc*2+i]
    return 100*fm/ff

def nx(iomat,i):
    im = 0
    ex = 0
    for j in range(0,nc):
        if j != i:
            ex = ex + iomat[i,j]+iomat[i,nc+j]+iomat[i,2*nc+j]
            im = im + iomat[j,i]+iomat[j,nc+i]+iomat[j,2*nc+i]
    return 100*(ex-im)/sum(iomat[4,:])
    #return 100*(ex-im)/iomat[4,i]

iomats=[]
for name in names:
    iomat2 = read_iomat(excel_path+'iomat-'+name+'.csv')
    iomats.append(iomat2)

with open(texpath + 'iomats-stats.tex','wb') as file:
    file.write('\\begin{table}[p]\n')
    file.write('\\renewcommand{\\arraystretch}{1.2}\n')
    file.write('\\begin{center}\n')
    file.write('\\caption{Measures of global production structure: benchmark and counterfactual input-output tables}\n')
    file.write('\\label{tab:iomats-stats}\n')
    #file.write('\\small\n')
    file.write('\\begin{tabular}{lcccc}\n')
    file.write('\\toprule\n')
    file.write('Region & ' + 
               '\\multicolumn{1}{p{3.5cm}}{\\centering World GDP share\\\\(percent)} & ' +
               '\\multicolumn{1}{p{3.5cm}}{\\centering Total trade \\\\(percent world GDP)} & ' +
               '\\multicolumn{1}{p{3.5cm}}{\\centering Intermediate trade \\\\(percent total trade)} & ' +
               '\\multicolumn{1}{p{3.5cm}}{\\centering Trade balance \\\\(percent world GDP)}\\\\\n')
    file.write('\\midrule\n')

    file.write('\\multicolumn{5}{l}{\\textit{(a) 1995 benchmark}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[0],i),
                    openness(iomats[0],i),
                    itrd(iomats[0],i),
                    nx(iomats[0],i)))
    file.write('\\midrule\n')

    file.write('\\multicolumn{5}{l}{\\textit{(b) 2011 benchmark}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[1],i),
                    openness(iomats[1],i),
                    itrd(iomats[1],i),
                    nx(iomats[1],i)))
    file.write('\\midrule\n')

    file.write('\\multicolumn{5}{l}{\\textit{(c) Counterfactual 1: change in relative sizes only}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[2],i),
                    openness(iomats[2],i),
                    itrd(iomats[2],i),
                    nx(iomats[2],i)))
    file.write('\\midrule\n')

    file.write('\\multicolumn{5}{l}{\\textit{(d) Counterfactual 2: change in openness only}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[3],i),
                    openness(iomats[3],i),
                    itrd(iomats[3],i),
                    nx(iomats[3],i)))
    file.write('\\midrule\n')

    file.write('\\multicolumn{5}{l}{\\textit{(e) Counterfactual 3: change in intermediate trade share only}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[4],i),
                    openness(iomats[4],i),
                    itrd(iomats[4],i),
                    nx(iomats[4],i)))
    file.write('\\midrule\n')

    file.write('\\multicolumn{5}{l}{\\textit{(f) Counterfactual 4: change in trade balances only}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[5],i),
                    openness(iomats[5],i),
                    itrd(iomats[5],i),
                    nx(iomats[5],i)))

    file.write('\\bottomrule\n')
    file.write('\\end{tabular}\n')
    file.write('\\normalsize\n')
    file.write('\\end{center}\n')
    file.write('\\end{table}\n')

# balanced trade version
with open(texpath + 'iomats-bal-stats.tex','wb') as file:
    file.write('\\begin{table}[p]\n')
    file.write('\\renewcommand{\\arraystretch}{1.2}\n')
    file.write('\\begin{center}\n')
    file.write('\\caption{Measures of global production structure in input-output tables}\n')
    file.write('\\label{tab:iomats-stats}\n')
    #file.write('\\small\n')
    file.write('\\begin{tabular}{lccc}\n')
    file.write('\\toprule\n')
    file.write('Region & ' + 
               '\\multicolumn{1}{p{3.5cm}}{\\centering World GDP share\\\\(percent)} & ' +
               '\\multicolumn{1}{p{3.5cm}}{\\centering Total trade \\\\(percent world GDP)} & ' +
               '\\multicolumn{1}{p{3.5cm}}{\\centering Intermediate trade \\\\(percent total trade)}\\\\\n')
    file.write('\\midrule\n')

    file.write('\\multicolumn{4}{l}{\\textit{(a) 1995 benchmark}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[6],i),
                    openness(iomats[6],i),
                    itrd(iomats[6],i)))
    file.write('\\midrule\n')

    file.write('\\multicolumn{4}{l}{\\textit{(b) 2011 benchmark}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[7],i),
                    openness(iomats[7],i),
                    itrd(iomats[7],i)))
    file.write('\\midrule\n')

    file.write('\\multicolumn{4}{l}{\\textit{(c) Counterfactual 1: change in size only}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[8],i),
                    openness(iomats[8],i),
                    itrd(iomats[8],i)))
    file.write('\\midrule\n')

    file.write('\\multicolumn{4}{l}{\\textit{(d) Counterfactual 2: change in trade openness only}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[9],i),
                    openness(iomats[9],i),
                    itrd(iomats[9],i)))
    file.write('\\midrule\n')

    file.write('\\multicolumn{4}{l}{\\textit{(e) Counterfactual 3: change in intermediate trade share only}}\\\\\n')
    for i in range(0,4):
        file.write(countries[i] + '& %0.2f & %0.2f & %0.2f\\\\\n' % 
                   (gdp_share(iomats[10],i),
                    openness(iomats[10],i),
                    itrd(iomats[10],i)))

    file.write('\\bottomrule\n')
    file.write('\\end{tabular}\n')
    file.write('\\normalsize\n')
    file.write('\\end{center}\n')
    file.write('\\end{table}\n')
