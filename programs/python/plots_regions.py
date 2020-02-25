##################################################################################
# This script merges WIOD data with portfolio diversification data, then creates
# the figures in the paper

##################################################################################
# Imports, constants, etc.

import regions
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.formula.api as sm

inpath1= 'results/'
inpath2= 'results/'
outpath= 'results/'
figpath= 'tabfig/'

LINEWIDTH = 3
COLORS = ['blue','red','green','purple','black']
#MARKERS = ['o','s','^','D',',']
MARKERS = [None,None,None,None,None,None]
#LINESTYLES = ['-','--','-.',':','-']
LINESTYLES = ['--','--','--','--','--','--']
FONTSIZE = 22
DASHES = [(None,None),(30,5),(10,3),(3,1),(1,0.5)]

mpl.rc('font',**{'family':'serif','serif':['Palatino'],'size':FONTSIZE})
mpl.rc('font',size=FONTSIZE)
mpl.rc('text', usetex=True)
mpl.rc('lines',linewidth=LINEWIDTH)

world = ['WLD']
maxyr=2011
years = range(1995,2012)

##################################################################################
# Function definitions

def tsplot(grouped,col,fname):
    # Creates a PDF of a time series plot with multiple lines from a Pandas GroupBy
    # object. Each line contains data in the specified column for a particular group.
    # Note this object must also have a column called year.

    legend=[]
    cnt=0
    for num,data in grouped:
        legend = legend + [regions.which_region_name(num)]
        plt.plot(data['year'],
                 data[col],
                 color=COLORS[cnt],
                 marker=MARKERS[cnt],
                 linestyle=LINESTYLES[cnt],
                 dashes=DASHES[cnt])
        cnt = cnt+1

    plt.xlim([1995,maxyr])
    plt.xticks([1995,1999,2003,2007,2011])
    plt.grid(False)
    plt.xlabel('')
    plt.tight_layout()
    plt.legend(legend,loc='best',prop={'size':16})
    plt.savefig(figpath+fname,bbox_inches='tight')
    plt.clf()
    plt.close()

##################################################################################
# Production structure plots

data = pd.read_csv(inpath1 + 'proddata_regions.txt')
data['region_num'] = data['region'].apply(lambda x: regions.which_region_num(x))
wld = data[data.region=='WLD'][['year','gdp']].rename(columns={'gdp':'wld_gdp'})
data = pd.merge(left=data,right=wld,how='left',on='year')
data['open'] = 100*data.gtrd/data.wld_gdp
data['tby'] = 100*data.tbe/data.wld_gdp
data['yshare'] = 100*data.gdp/data.wld_gdp
data['ishare'] = 100*data.itrd/data.gtrd
data['ishare_im'] = 100*data.imi/data.im
data['ishare_ex'] = 100*data.exi/data.ex
data['mshare_go'] = 100*(data.go-data.gdp)/data.go
data = data[data.region != 'WLD']

grouped = data.groupby('region_num')
tsplot(grouped,'open','ts_openness.pdf')
tsplot(grouped,'ishare','ts_itrade.pdf')
tsplot(grouped,'tby','ts_trade_balances.pdf')
tsplot(grouped,'yshare','ts_gdp_shares.pdf')

##################################################################################
# Portfolios plots

data = pd.read_csv(inpath2 + 'portfolios.csv')
data = data[data.year>=1995]
data = data[pd.notnull(data.country)]
data['region'] = data['country'].apply(lambda x: regions.which_region(x))

# group by region and year
data['divprod1'] = data['weight1']*data['div']/100
data['divprod2'] = data['weight2']*data['div']/100

data['far_prod1'] = data['weight1']*data['fa_risky_frac']/100
data['far_prod2'] = data['weight2']*data['fa_risky_frac']/100
data['flr_prod1'] = data['weight1']*data['fl_risky_frac']/100
data['flr_prod2'] = data['weight2']*data['fl_risky_frac']/100
data['fas_prod1'] = data['weight1']*data['fa_safe_frac']/100
data['fas_prod2'] = data['weight2']*data['fa_safe_frac']/100
data['fls_prod1'] = data['weight1']*data['fl_safe_frac']/100
data['fls_prod2'] = data['weight2']*data['fl_safe_frac']/100
data['fdia_prod1'] = data['weight1']*data['fdia_frac']/100
data['fdia_prod2'] = data['weight2']*data['fdia_frac']/100
data['fdil_prod1'] = data['weight1']*data['fdil_frac']/100
data['fdil_prod2'] = data['weight2']*data['fdil_frac']/100
data['rfa_prod1'] = data['weight1']*data['risky_frac_assets']/100
data['rfa_prod2'] = data['weight2']*data['risky_frac_assets']/100
data['sfl_prod1'] = data['weight1']*data['safe_frac_liabilities']/100
data['sfl_prod2'] = data['weight2']*data['safe_frac_liabilities']/100

grouped = data.groupby(['region','year'])

# aggregate 
cols = ['weight1',
        'weight2',
        'divprod1',
        'divprod2',
        'far_prod1',
        'far_prod2',
        'flr_prod1',
        'flr_prod2',
        'fas_prod1',
        'fas_prod2',
        'fls_prod1',
        'fls_prod2',
        'fdia_prod1',
        'fdia_prod2',
        'fdil_prod1',
        'fdil_prod2',
        'rfa_prod1',
        'rfa_prod2',
        'sfl_prod1',
        'sfl_prod2']

summed = grouped[cols].sum()
summed = summed.reset_index()

summed['div_mean1'] = 100*summed['divprod1']/summed['weight1']
summed['div_mean2'] = 100*summed['divprod2']/summed['weight2']

summed['far_mean1'] = 100*summed['far_prod1']/summed['weight1']
summed['far_mean2'] = 100*summed['far_prod2']/summed['weight2']
summed['flr_mean1'] = 100*summed['flr_prod1']/summed['weight1']
summed['flr_mean2'] = 100*summed['flr_prod2']/summed['weight2']
summed['fas_mean1'] = 100*summed['fas_prod1']/summed['weight1']
summed['fas_mean2'] = 100*summed['fas_prod2']/summed['weight2']
summed['fls_mean1'] = 100*summed['fls_prod1']/summed['weight1']
summed['fls_mean2'] = 100*summed['fls_prod2']/summed['weight2']
summed['fdia_mean1'] = 100*summed['fdia_prod1']/summed['weight1']
summed['fdia_mean2'] = 100*summed['fdia_prod2']/summed['weight2']
summed['fdil_mean1'] = 100*summed['fdil_prod1']/summed['weight1']
summed['fdil_mean2'] = 100*summed['fdil_prod2']/summed['weight2']
summed['rfa_mean1'] = 100*summed['rfa_prod1']/summed['weight1']
summed['rfa_mean2'] = 100*summed['rfa_prod2']/summed['weight2']
summed['sfl_mean1'] = 100*summed['sfl_prod1']/summed['weight1']
summed['sfl_mean2'] = 100*summed['sfl_prod2']/summed['weight2']

medians = grouped[['div',
                   'fa_risky_frac',
                   'fl_risky_frac',
                   'fa_safe_frac',
                   'fl_safe_frac',
                   'fdia_frac',
                   'fdil_frac',
                   'risky_frac_assets',
                   'safe_frac_liabilities']].agg(np.median)
medians = medians.reset_index()
medians = medians.rename(columns={'div':'div_median',
                                  'fa_risky_frac':'far_median',
                                  'fl_risky_frac':'flr_median',
                                  'fa_safe_frac':'fas_median',
                                  'fl_safe_frac':'fls_median',
                                  'fdia_frac':'fdia_median',
                                  'fdil_frac':'fdil_median',
                                  'risky_frac_assets':'rfa_median',
                                  'safe_frac_liabilities':'sfl_median'})
data2 = pd.merge(left=summed,right=medians,how='left',on=['region','year'])
data2['region_num'] = data2['region'].apply(lambda x: regions.which_region_num(x))
data2 = data2.sort(['region_num','year']).reset_index(drop=True)

grouped = data2.groupby('region_num')
tsplot(grouped,'div_mean1','ts_div_mean1.pdf')
tsplot(grouped,'div_mean2','ts_div_mean2.pdf')
tsplot(grouped,'div_median','ts_div_median.pdf')

all_plots=0
if(all_plots==1):
    tsplot(grouped,'far_mean1','ts_far_mean1.pdf')
    tsplot(grouped,'far_mean2','ts_far_mean2.pdf')
    tsplot(grouped,'far_median','ts_far_median.pdf')
    
    tsplot(grouped,'flr_mean1','ts_flr_mean1.pdf')
    tsplot(grouped,'flr_mean2','ts_flr_mean2.pdf')
    tsplot(grouped,'flr_median','ts_flr_median.pdf')

    tsplot(grouped,'fas_mean1','ts_fas_mean1.pdf')
    tsplot(grouped,'fas_mean2','ts_fas_mean2.pdf')
    tsplot(grouped,'fas_median','ts_fas_median.pdf')

    tsplot(grouped,'fls_mean1','ts_fls_mean1.pdf')
    tsplot(grouped,'fls_mean2','ts_fls_mean2.pdf')
    tsplot(grouped,'fls_median','ts_fls_median.pdf')

    tsplot(grouped,'fdia_mean1','ts_fdia_mean1.pdf')
    tsplot(grouped,'fdia_mean2','ts_fdia_mean2.pdf')
    tsplot(grouped,'fdia_median','ts_fdia_median.pdf')

    tsplot(grouped,'fdil_mean1','ts_fdil_mean1.pdf')
    tsplot(grouped,'fdil_mean2','ts_fdil_mean2.pdf')
    tsplot(grouped,'fdil_median','ts_fdil_median.pdf')

    tsplot(grouped,'rfa_mean1','ts_rfa_mean1.pdf')
    tsplot(grouped,'rfa_mean2','ts_rfa_mean2.pdf')
    tsplot(grouped,'rfa_median','ts_rfa_median.pdf')

    tsplot(grouped,'sfl_mean1','ts_sfl_mean1.pdf')
    tsplot(grouped,'sfl_mean2','ts_sfl_mean2.pdf')
    tsplot(grouped,'sfl_median','ts_sfl_median.pdf')
 
# last make a datafile for matlab comparisons
first = data2[data2.year==1995].reset_index(drop=True)
last = data2[data2.year==2011].reset_index(drop=True)
delta = last.copy()
for c in first.columns:
    if c!='region_num' and c!='region':
        delta[c] = delta[c] - first[c]

cols = ['div_median','div_mean1','div_mean2']
with open(outpath +'div_for_matlab.csv','wb') as file:
    file.write('region_num,med_0,mean1_0,mean2_0,med_1,mean1_1,mean2_1,med_d,mean1_d,mean2_d\n')
    for i in range(regions.nr):
        file.write('%s,' % first['region_num'][i])
        for c in cols:
            file.write('%f,' % first[c][i])
        for c in cols:
            file.write('%f,' % last[c][i])
        for c in cols:
            file.write('%f,' % delta[c][i])
        file.write('\n')
