import numpy as np
import pandas as pd
import scikits.statsmodels.tsa.api as sm
import statsmodels
import regions

# load the data
inpath = '/media/data/datasets/pwt/'
outpath = '/home/joe/my-research/prodchains/programs/matlab/portfolios/wiod-data/adv-eme/'

pwt = pd.read_stata(inpath+'pwt80.dta')
pwt = pwt.sort(['countrycode','year'])

# filter tfp by country (using rtfpna), then run AR1 on each country separately
lam=500
def hpcycle(x):
    cycle,trend = sm.filters.hpfilter(x,lamb=lam)
    return cycle,trend

pwt2 = pwt[pd.notnull(pwt.rtfpna)]
ccodes = pwt2.countrycode.unique()

pwt2['ltfp'] = np.log(pwt2.rtfpna)
pwt2['ltfp-cycle'] = 0.0
pwt2['ltfp-trend'] = 0.0
res = pd.DataFrame(columns=['countrycode','region','nobs','rho_z','sigma_z'],index=range(len(ccodes)))
res['rho_z']=0.0
res['sigma_z']=0.0
cnt=0
for c,d in pwt2.groupby('countrycode'):
    index = d.index
    cc,tt=hpcycle(d.ltfp)
    pwt2['ltfp-cycle'][index]=cc
    pwt2['ltfp-trend'][index]=tt

    res['countrycode'][cnt]=c
    res['region'][cnt]=regions.which_region(c)
    res['nobs'][cnt]=len(cc)
    model = sm.AR(cc)
    results=model.fit(1)
    res['rho_z'][cnt]=results.params[1]
    res['sigma_z'][cnt]=np.std(results.resid)    
    cnt=cnt+1

grouped = res.groupby('region')
process = grouped[['rho_z','sigma_z']].agg(np.median)

# aggregate across regions first, then run VAR(1) to get joint tfp process
mask = np.logical_and(pd.notnull(pwt.rgdpo),pd.notnull(pwt.emp))
mask = np.logical_and(mask,pd.notnull(pwt.avh))
mask = np.logical_and(mask,pd.notnull(pwt.rkna))
mask = np.logical_and(mask,pd.notnull(pwt.ck))
mask = np.logical_and(mask,pd.notnull(pwt.cgdpo))
mask = np.logical_and(mask,pd.notnull(pwt.hc))
pwt2=pwt[mask]

pwt2['y'] = pwt2.rgdpo
pwt2['k'] = pwt2.rgdpo*(pwt2.ck/pwt2.cgdpo)
pwt2['n'] = pwt2.emp*pwt2.avh*pwt.hc

alpha=0.32
pwt2['region'] = pwt2.countrycode.apply(lambda x: regions.which_region(x))
pwt2['region_num'] = pwt2.region.apply(lambda x: regions.which_region_num(x))
grouped = pwt2.groupby(['region_num','year'])
agged = grouped[['k','n','y']].sum()
agged=agged.reset_index()
agged['ltfp'] = np.log(agged.y) - alpha*np.log(agged.k) - (1-alpha)*np.log(agged.n)

agged['ltfp-cycle']=0.0
agged['ltfp-trend']=0.0
agged=agged.sort(['region_num','year'])
for r,d in agged.groupby('region_num'):
    index=d.index
    cc,tt=hpcycle(d.ltfp)
    agged['ltfp-cycle'][index]=cc
    agged['ltfp-trend'][index]=tt

# correlated shocks
pivoted = agged[['region_num','year','ltfp-cycle']].pivot('year','region_num')
model=sm.VAR(pivoted.values)
results=model.fit(1)
N=results.coefs
SIGMA=np.cov(np.transpose(results.resid))
np.savetxt(outpath+'N.txt',N[0],fmt='%0.15f',delimiter=' ')
np.savetxt(outpath+'SIGMA.txt',SIGMA,fmt='%0.15f',delimiter=' ')

# uncorrelated shocks
N2 = np.zeros((4,4))
SIG2 = np.zeros((4,4))
for i in range(4):
    model = sm.AR(pivoted.values[:,i])
    results=model.fit(1)
    N2[i,i]=results.params[1]
    SIG2[i,i]=(np.std(results.resid))**2.0
np.savetxt(outpath+'N2.txt',N2,fmt='%0.15f',delimiter=' ')
np.savetxt(outpath+'SIG2.txt',SIG2,fmt='%0.15f',delimiter=' ')
