##########################################################################################
# This script calculates a measure of portfolio diversification following the methodology
# of Heathcote and Perri (JPE 2014).
#
# Data sources:
# gross foreign assets (FA): Lane and Milesi-Feretti
# gross foreign liabilities (FL): Lane and Miesli-Feretti
# capital stock (K): constructed using perpetual inventory method from PWT8.0
##########################################################################################

##########################################################################################
# Imports, constants, etc.
##########################################################################################
import numpy as np
import pandas as pd
import regions

# paths
inpath = '/home/joe/Datasets/'
outpath = 'results/'

##########################################################################################
# Main script body
##########################################################################################

# Load Lane and Milesi-Feretti data and calculate FA and FL as fractions of GDP
# Also calculate risky and safe portions of FA and FL

lmf = pd.read_csv(inpath+'lane_milesi-ferretti/ewn19702011.csv',sep=',')
lmf = lmf.rename(columns=lambda x: x.strip())
lmf = lmf.rename(columns={'Year':'year'})

assets = ['Portfolio equity assets (stock)',
          'FDI assets (stock)',
          'Debt assets (stock)',
          'financial derivatives (assets)',
          'FX Reserves minus gold']
assets_risky = ['Portfolio equity assets (stock)',
                'FDI assets (stock)',
                'financial derivatives (assets)']
assets_safe = ['Debt assets (stock)','FX Reserves minus gold']

lmf['fa'] = 0
lmf['fa'] = lmf['fa'].astype(float)
for a in assets:
    lmf['fa'] = lmf['fa'] + lmf[a]

lmf['fa_risky'] = 0
lmf['fa_risky'] = lmf['fa_risky'].astype(float)
for a in assets_risky:
    lmf['fa_risky'] = lmf['fa_risky'] + lmf[a]

lmf['fa_safe'] = 0
lmf['fa_safe'] = lmf['fa_safe'].astype(float)
for a in assets_safe:
    lmf['fa_safe'] = lmf['fa_safe'] + lmf[a]

lmf['fdia'] = lmf['FDI assets (stock)'].astype(float)

liabilities = ['Portfolio equity liabilities (stock)',
               'FDI liabilities (stock)',
               'Debt liabilities (stock)',
               'financial derivatives (liab)']
liabilities_risky = ['Portfolio equity liabilities (stock)',
                      'FDI liabilities (stock)',
                      'financial derivatives (liab)']
liabilities_safe = ['Debt liabilities (stock)']


lmf['fl'] = 0
lmf['fl'] = lmf['fl'].astype(float)
for l in liabilities:
    lmf['fl'] = lmf['fl'] + lmf[l]

lmf['fl_risky'] = 0
lmf['fl_risky'] = lmf['fl_risky'].astype(float)
for l in liabilities_risky:
    lmf['fl_risky'] = lmf['fl_risky'] + lmf[l]

lmf['fl_safe'] = 0
lmf['fl_safe'] = lmf['fl_safe'].astype(float)
for l in liabilities_safe:
    lmf['fl_safe'] = lmf['fl_safe'] + lmf[l]

lmf['fdil'] = lmf['FDI liabilities (stock)'].astype(float)

lmf['fay'] = lmf['fa']/lmf['GDP (US$)']
lmf['fly'] = lmf['fl']/lmf['GDP (US$)']
lmf['fay_risky'] = lmf['fa_risky']/lmf['GDP (US$)']
lmf['fly_risky'] = lmf['fl_risky']/lmf['GDP (US$)']
lmf['fay_safe'] = lmf['fa_safe']/lmf['GDP (US$)']
lmf['fly_safe'] = lmf['fl_safe']/lmf['GDP (US$)']
lmf['fdiay'] = lmf['fdia']/lmf['GDP (US$)']
lmf['fdily'] = lmf['fdil']/lmf['GDP (US$)']

####################################################################

# Load PWT data
pwt = pd.read_stata(inpath+'pwt/pwt80.dta')
pwt['ky1'] = pwt['ck']/pwt['cgdpo']

# Concordance between country codes and ifs codes
clmf = lmf[['Country Name','IFS id']].drop_duplicates()
cpwt = pwt[['country','countrycode']].drop_duplicates()
clmf.to_csv(outpath+'clmf.csv',sep=',',encoding='utf-8')
cpwt.to_csv(outpath+'cpwt.csv',sep=',',encoding='utf-8')
# Note: at this stage we manually construct the concordance and store it in countries.csv
countries = pd.read_csv(outpath + 'countries.csv',sep='\t')
merged1 = pd.merge(left=lmf,right=countries,how='left',on='IFS id')
merged2 = pd.merge(left=merged1,right=pwt[['countrycode','year','ky1']],how='left',on=['countrycode','year'])

# calculate diversification measure
merged2['div'] = 100*(merged2['fay']+merged2['fly'])/(2*(merged2['ky1']+merged2['fay']-merged2['fly']))
merged2['fa_risky_frac'] = 100*(merged2['fay_risky'])/(merged2['ky1']+merged2['fay']-merged2['fly'])
merged2['fl_risky_frac'] = 100*(merged2['fly_risky'])/(merged2['ky1']+merged2['fay']-merged2['fly'])
merged2['fa_safe_frac'] = 100*(merged2['fay_safe'])/(merged2['ky1']+merged2['fay']-merged2['fly'])
merged2['fl_safe_frac'] = 100*(merged2['fly_safe'])/(merged2['ky1']+merged2['fay']-merged2['fly'])
merged2['fdia_frac'] = 100*merged2['fdiay']/(merged2['ky1']+merged2['fay']-merged2['fly'])
merged2['fdil_frac'] = 100*merged2['fdily']/(merged2['ky1']+merged2['fay']-merged2['fly'])
merged2['risky_frac_assets'] = 100*merged2['fay_risky']/(merged2['fay_risky']+merged2['fay_safe'])
merged2['safe_frac_liabilities'] = 100*merged2['fly_safe']/(merged2['fly_risky']+merged2['fly_safe'])
#merged2['weight1'] = merged2['fa']-merged2['fl']+merged2['ky1']*merged2['GDP (US$)']
merged2['weight1'] = merged2['ky1']*merged2['GDP (US$)']
merged2['weight2'] = merged2['GDP (US$)']

# countries we do not want to include in the DIV data because they are craxy outliers
mask = merged2['countrycode'].isin([
    'PAN', # panama
    'SYR', # syria
    'VNM', # vietnam
    'YEM', # yemen
    #'CYP', # cyprus
    #'GBR', # great britain
    'AGO', # angola
    'ATG', # antigua
    'BHR', # bahrain
    'BLZ', # belize
    #'BEN', # benin
    'BIH', # !! incomplete data, 1998 onward only
    'HRV', # !! incomplete data, 1996 onward only
    'TCD', # chad
    'COG', # congo
    'CIV', # cote d'ivoire
    #'CHN', # china; capital controls
    'DMA', # dominica
    #'EGY', # egypt
    'GNQ', # equatorial guinea
    'ETH', # ethiopia
    #'GAB', # gabon
    #'GMB', # gambia
    'GRD', # grenada
    'ISL', # iceland
    'IRL', # ireland
    'LAO', # laos
    'LBN', # lebanon
    'LBR', # liberia
    'LUX', # luxembourg
    'MDG', # madagascar
    'MLI', # mali
    'MOZ', # mozambique
    'NGA', # nigeria
    'STP', # sao tome
    'SLE', # sierra leone
    'LCA', # st. lucia
    'VCT', # st vincent
    'SDN', # sudan
    'TTO', # trinidad & tobago
    'ZMB', # zambia
    'ZWE' # zimbabwe
])

merged2.loc[mask,'weight1'] = 0.0
merged2.loc[mask,'weight2'] = 0.0

# write to file
merged2=merged2.rename(columns={'countrycode':'country'})
merged2[['country',
         'year',
         'div',
         'fa_risky_frac',
         'fl_risky_frac',
         'fa_safe_frac',
         'fl_safe_frac',
         'fdia_frac',
         'fdil_frac',
         'risky_frac_assets',
         'safe_frac_liabilities',
         'weight1',
         'weight2']].to_csv(outpath+'portfolios.csv',sep=',')
