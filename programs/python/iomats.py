##################################################################################
# This script does the following:
# - Constructs the raw, unbalanced IO matrices from processed WIOD data files
# - Writes this matrices to csv files (which I balance in Excel manually)
# - Reads in the balanced benchmark IO matrices, as well as the alternatives creates in Excel
# - Writes all of these balanced matrices to files used by MATLAB to calibrate parameters
# - Writes latex tables for each matrix

##################################################################################
# Imports, constants, etc.

import pandas as pd
import numpy as np
import itertools

ns=1
nc=4

years=range(1995,2012)
excel_path='../../excel/'
cal_path = 'calibration_data/'
texpath = 'tabfig/'
inpath = 'wiod_data/'
countries = ['USA','ADV','EME','ROW']
outpath= 'wiod_data/'

##################################################################################
# Function definitions

def read_wiod_data(year):
    # Reads the WIOD data for a given year and returns three Pandas DataFrames:
    # va: value added
    # inin: intermediate inputs
    # fin: final demand
    va = pd.read_csv(inpath + 'va' + str(year) +'.txt',sep=' ')
    inin = pd.read_csv(inpath + 'inin' + str(year) +'.txt',sep=' ')
    fin = pd.read_csv(inpath + 'fin' + str(year) +'.txt',sep=' ')
    va=va[pd.notnull(va['i'])]
    inin=inin[pd.notnull(inin['i'])]
    fin=fin[pd.notnull(fin['i'])]
    return va, inin, fin

def read_iomat(fname):
    # Reads an IO matrix from text file and returns a numpy matrix
    iomat = np.genfromtxt(fname=fname,dtype='float',delimiter=',',names=None)
    rowsums=np.zeros(nc*ns+1)
    colsums=np.zeros(nc*ns+nc*2)

    for row in range(0,nc*ns + 1):
        rowsums[row] = np.sum(iomat[row,:])

    for col in range(0,nc*ns + nc*2):
        colsums[col] = np.sum(iomat[:,col])

    return iomat, rowsums, colsums

def construct_iomat(va,inin,fin):
    # Given three Pandas DataFrames that contain value added (va), intermediate inputs (inin),
    # and final demand (fin), constructs a world input-output matrix. Returns three Numpy arrays:
    # iomat: IO matrix
    # rowsums: sums of rows in iomat
    # colsums: sums of columns in iomat
    va=va.sort(['i','s'])
    fin=fin.sort(['i','j','r'])
    inin=inin.sort(['i','s','j','r'])

    rowsums = np.zeros( nc*ns +1 )
    colsums = np.zeros( nc*ns + nc*2 )

    M = inin.pivot_table(values='m', rows=['j','r'], cols=['i','s']).values
    V = va['va'].values.reshape((1,nc*ns))
    F = fin.pivot_table(values=['c','x'], rows=['j','r'], cols='i').values
    V = np.hstack((V,np.zeros((1,nc*2))))
    iomat=np.vstack( ( np.hstack((M,F)) , V ) )

    for row in range(0,nc*ns + 1):
        rowsums[row] = np.sum(iomat[row,:])

    for col in range(0,nc*ns + nc*2):
        colsums[col] = np.sum(iomat[:,col])

    return iomat, rowsums, colsums

def deconstruct_iomat(iomat,rowsums,colsums):
    # Given a world input-output matrix (iomat), rowsums, and colsums, deconstructs it
    # using the opposite procedure as construct_iomat. Returns three Pandas DataFrames:
    # va: value added
    # inin: intermediate inputs
    # fin: final demand
    M=iomat[0:(nc*ns),0:(nc*ns)]
    V=iomat[-1,0:(nc*ns)]
    Fc=iomat[0:(nc*ns),(nc*ns):+((nc*ns)+nc)]
    Fx=iomat[0:(nc*ns),((nc*ns)+nc):]
    x=list(itertools.product(range(0,nc),range(0,ns)))

    va=pd.DataFrame(index=range(0,nc*ns),columns=['i','s','va','go'])
    for t in range(0,nc*ns):
        va['i'][t]=x[t][0]+1
        va['s'][t]=x[t][1]+1
    va['va'] = V
    va['go'] = colsums[0:(nc*ns)]
    
    inin=pd.DataFrame(index=range(0,nc*ns),columns=['j','r']+['m_'+str(t) for t in range(0,(nc*ns))])
    for t in range(0,nc*ns):
        inin.loc[t,'j']=x[t][0]+1
        inin.loc[t,'r']=x[t][1]+1
        inin['m_'+str(t)] = M[:,t]
    inin=pd.melt(inin,id_vars=['j','r'])
    inin=inin.rename(columns={'value':'m'})
    y=dict(zip(['m_'+str(t) for t in range(0,(nc*ns))],x))
    inin['i']=inin['variable'].apply(lambda x: y[x][0]+1)
    inin['s']=inin['variable'].apply(lambda x: y[x][1]+1)    
    inin=inin[['i','s','j','r','m']]

    finc=pd.DataFrame(index=range(0,nc*ns),columns=['j','r']+['c_'+str(t) for t in range(0,nc)])
    finx=pd.DataFrame(index=range(0,nc*ns),columns=['j','r']+['x_'+str(t) for t in range(0,nc)])
    for t in range(0,nc*ns):
        finc['j'][t]=x[t][0]+1
        finx['j'][t]=x[t][0]+1
        finc['r'][t]=x[t][1]+1
        finx['r'][t]=x[t][1]+1
    for t in range(0,nc):
        finc['c_'+str(t)] = Fc[:,t]
        finx['x_'+str(t)] = Fx[:,t]
    finc=pd.melt(finc,id_vars=['j','r'])
    finx=pd.melt(finx,id_vars=['j','r'])
    finc=finc.rename(columns={'value':'c'})
    finx=finx.rename(columns={'value':'x'})
    yc=dict(zip(['c_'+str(t) for t in range(0,nc)],range(0,nc)))
    yx=dict(zip(['x_'+str(t) for t in range(0,nc)],range(0,nc)))
    finc['i']=finc['variable'].apply(lambda x: yc[x]+1)
    finx['i']=finx['variable'].apply(lambda x: yx[x]+1)

    finc=finc[['i','j','r','c']]
    finx=finx[['i','j','r','x']]
    fin=pd.merge(left=finc,right=finx,how='left',on=['i','j','r'])

    return va, inin, fin

def write_iomat_csv(iomat,fname):
    # Write world IO matrix (iomat) to csv file called filename
    np.savetxt(fname=outpath+fname,X=iomat,fmt='%0.15f',delimiter=',')
    
def write_iomat_latex(iomat,rowsums,colsums,caption,label,fname):
    # Given a world IO matrix (iomat), rowsums, colsums, creates a latex file
    # in location filename that contains a table with given caption and label.

    iomat2=100*iomat[:,:]/rowsums[-1]
    rowsums2=100*rowsums/rowsums[-1]
    colsums2=100*colsums/rowsums[-1]
    M=iomat2[0:(nc*ns),0:(nc*ns)]
    V=iomat2[-1,0:(nc*ns)]
    Fc=iomat2[0:(nc*ns),(nc*ns):+((nc*ns)+nc)]
    Fx=iomat2[0:(nc*ns),((nc*ns)+nc):]

    with open(texpath + fname + '.tex','wb') as file:
        file.write('\\begin{table}[p]\n')
        file.write('\\begin{center}\n')
        file.write('\\caption{'+caption+'}\n')
        file.write('\\label{tab:'+label+'}\n')
        file.write('\\small\n')

        file.write('\\begin{tabular}{c')
        for i in range(0,nc):
            file.write('cc')
        file.write('c}\n')
        file.write('\\toprule\n')
            
        file.write('& \\multicolumn{'+str(nc)+'}{c}{Intermediate inputs}')
        file.write('& \\multicolumn{'+str(nc)+'}{c}{Final demand} & \\\\\n')
        file.write('\\cmidrule(rl){2-'+str(nc+1)+'}\\cmidrule(rl){'+str(nc+2)+'-'+str(nc+1+nc)+'}\n')
        for c in countries:
            file.write(' &'+c)
        for c in countries:
            file.write(' &'+c)
        file.write('& GO\\\\\n')
        file.write('\\midrule\n')

        for i in range(0,nc):
            file.write(countries[i])
            for j in range(0,nc):
                file.write('& %0.2f' % M[i,j])
            for j in range(0,nc):
                file.write('& %0.2f' % (Fc[i,j]+Fx[i,j]))
            file.write('& %0.2f '% rowsums2[i])
            file.write('\\\\\n')
            
        file.write('\\midrule\n')
        file.write('VA')
        for i in range(0,nc):
            file.write('& %0.2f' % V[i])
        for i in range(0,nc):
            file.write('& 0.00')
        file.write('& %0.2f' % rowsums2[-1])
        file.write('\\\\\n')

        file.write('\\midrule\n')
        file.write('GO')
        for i in range(0,nc):
            file.write('& %0.2f' % colsums2[i])
        for i in range(0,nc):
            file.write('& 0.00')
        file.write('& %0.2f' % sum(colsums2[0:nc]))
        file.write('\\\\\n')
            
        file.write('\\bottomrule\n')
        file.write('\\end{tabular}\n')
        file.write('\\normalsize\n')
        file.write('\\end{center}\n')
        file.write('\\end{table}\n')

def coeffs(iomat):
    # Given world IO matrix (iomat), calculates IO coefficients and returs them in A
    A=np.zeros(iomat.shape)
    for col in range(0,A.shape[1]):
        A[:,col] = iomat[:,col]/np.sum(iomat[:,col])
    return A

def ras(iomat0,rowsums1,colsums1):
    # Given an initial IO matrix (iomat), and desired rowsums (rowsums1) and colsums (colsums1),
    # performs the RAS balancing procedure. Returns a new IO matrix (iomat) that is consistent
    # with desired row- and colsums.
    A0 = coeffs(iomat0)
    iomat = np.dot(A0,np.diag(colsums1))

    go=True
    iter=0
    maxit=1000
    tol=1.0e-6

    while go:
        iter=iter+1
        rowsums = np.sum(iomat,axis=1)
        r = np.divide(rowsums1,rowsums)
        iomat = np.dot(np.diag(r),iomat)
        colsums = np.sum(iomat,axis=0)
        s = np.divide(colsums1,colsums)
        iomat = np.dot(iomat,np.diag(s))
        colsums = np.sum(iomat,axis=0)
        rowsums = np.sum(iomat,axis=1)

        norm1 = max(np.divide(abs(rowsums-rowsums1),rowsums1))
        norm2 = max(np.divide(abs(colsums-colsums1),colsums1))
        if((norm1 <tol and norm2 <tol) or iter == maxit):
            go=False

        if iter==maxit:
            print 'RAS iteration did not converge!'
            print 'iter = ', iter, ' diff = ', max(norm1,norm2)

    return iomat

##################################################################################
# Main script body

# first load the raw WIOD data, create the (slightly) resulting unbalanced IO matrices,
# and write them to Excel so we can balance them and create the alternatives

load_wiod = False

if load_wiod==True:
    for year in years:
        va, inin, fin = read_wiod_data(year)
        iomat, rowsums, colsums = construct_iomat(va, inin, fin)
        write_iomat_csv(iomat,'iomat-raw-'+str(year)+'.csv')

# load the benchmark and alternative tables, make latex tables, and create the MATLAB calibration files
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

for name,title in zip(names,titles):
    
    # load the matrix
    iomat, rowsums, colsums = read_iomat(excel_path+'iomat-'+name+'.csv')

    # write the latex file
    write_iomat_latex(iomat,rowsums,colsums,title,'iomat-'+name,'iomat-'+name)

    # create calibration text files for matlab program to use
    va, inin, fin = deconstruct_iomat(iomat,rowsums,colsums)

    va.to_csv(cal_path + 'va-' + name + '.txt',
                 sep=' ',columns=['i','s','go','va'],index=False)
    inin.to_csv(cal_path + 'inin-'+ name + '.txt',
                   sep=' ',columns=['i','s','j','r','m'],index=False)
    fin.to_csv(cal_path + 'fin-' + name + '.txt',
                  sep=' ',columns=['i','j','r','c','x'],index=False)
    with open(cal_path + 'va-' + name + '.txt', 'a') as myfile:
        myfile.write('\n')
    with open(cal_path + 'inin-' + name + '.txt', 'a') as myfile:
        myfile.write('\n')
    with open(cal_path + 'fin-' + name + '.txt', 'a') as myfile:
        myfile.write('\n')
