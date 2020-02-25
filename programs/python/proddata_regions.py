##########################################################################################
# This script uses WIOD data to calculate various measures of vertical specialization,
# production fragmentation, etc.
#
# The resulting datafile contains the following variables:
# country: obvious
# year: obvious
# gdp: obvious
# go: gross output
# exy: gross exports as fraction of GDP
# imy: gross imports as fraction of GDP
# exf_y: gross final exports as fraction of GDP
# exi_y: gross intermediate exports as fraction of GDP
# imf_y: gross final imports as fraction of GDP
# imi_y: gross intermediate imports as fraction of GDP
# vxy: value added exports as fraction of GDP (Johnson and Noguera, 2012)
# tbe: gross trade balance/GDP
# tbv: value added trade balance/GDP
# dvx: domestic value added (KWW, 2014)
# dcx: domestic content (KWW, 2014)
# vsx: foreign value added (KWW, 2014; similar to Hummels et al., 2001)
# gtrdy: gross trade as fraction of GDP
# itrdy: gross intermediate trade as fraction of GDP
# ftrdy: gross final trade as fraction of GDP
# vtrdy: value added trade as fraction of GDP
# ishare: intermediate trade as fraction of total trade
##########################################################################################

##########################################################################################
# Imports, constants, etc.
##########################################################################################
import pickle
import pandas as pd
import numpy as np
import regions

years=range(1995,2012)
inpath= '/media/data_/datasets/WIOD/pik/'
outpath= 'results/'

iosys=[]

nr = 4

##########################################################################################
# Main script body
##########################################################################################

regions_ = regions.regions_dict.keys()
regions_.append('WLD')
go=np.zeros((len(years),nr+1))
gdp=np.zeros((len(years),nr+1))
ex=np.zeros((len(years),nr+1))
im=np.zeros((len(years),nr+1))
exf=np.zeros((len(years),nr+1))
exi=np.zeros((len(years),nr+1))
imf=np.zeros((len(years),nr+1))
imi=np.zeros((len(years),nr+1))
tbe=np.zeros((len(years),nr+1))
gtrd=np.zeros((len(years),nr+1))
itrd=np.zeros((len(years),nr+1))
ftrd=np.zeros((len(years),nr+1))
ishare=np.zeros((len(years),nr+1))

# for each year...
for t in range(len(years)):

    year=years[t]
    with open(inpath+'wiot_iosys_'+str(year)+'.pik','rb') as file:
        iosys=pickle.load(file)

    Y=iosys['Y']
    VAR=iosys['VAR']
    E=iosys['E']
    E_I=iosys['E_I']
    E_J=iosys['E_J']
    C=iosys['C']
    
    sumgo = 0
    sumgdp = 0
    sumex = 0
    sumim = 0
    sumex_i = 0
    sumex_f = 0
    sumim_i = 0
    sumim_f = 0
    sumgtrd = 0
    sumitrd = 0
    sumftrd = 0

    # and each region...
    for i in range(0,nr):
        
        r = regions_[i]

        # get the data where region val is the source and other regions are the dest
        masksrc = np.array([regions.which_region(e) == r for e in E_I],dtype=bool)
        maskdst = np.array([regions.which_region(e) != r for e in E_J],dtype=bool)

        ex2 = E[masksrc,:][:,maskdst].sum()
        ex_f = C[masksrc,:][:,maskdst].sum()
        ex_i = ex2 - ex_f
        gdp2 = np.multiply(VAR[masksrc],Y[masksrc]).sum()
        go2 = Y[masksrc].sum()
        
        sumex = sumex+ex2
        sumex_f = sumex_f + ex_f
        sumex_i = sumex_i + ex_i
        sumgdp = sumgdp + gdp2
        sumgo = sumgo + go2

        # now get the data where region val in the dest and other regions are the source
        masksrc = np.array([regions.which_region(e) != r for e in E_I],dtype=bool)
        maskdst = np.array([regions.which_region(e) == r for e in E_J],dtype=bool)

        im2 = E[masksrc,:][:,maskdst].sum()
        im_f = C[masksrc,:][:,maskdst].sum()
        im_i = im2 - im_f

        sumgtrd = sumgtrd + ex2+im2
        sumitrd = sumitrd + ex_i + im_i
        sumftrd = sumftrd + ex_f + im_f
        
        go[t,i] = go2
        gdp[t,i] = gdp2
        ex[t,i] = ex2
        im[t,i] = im2
        exi[t,i] = ex_i
        exf[t,i] = ex_f
        imi[t,i] = im_i
        imf[t,i] = im_f
        tbe[t,i] = ((ex2-im2))
        gtrd[t,i] = ((ex2+im2))
        itrd[t,i] = ((ex_i+im_i))
        ftrd[t,i] = ((ex_f+im_f))
        ishare[t,i] = ((ex_i+im_i)/(ex2+im2))

        
    go[t,-1] = sumgo
    gdp[t,-1] = sumgdp
    ex[t,-1] = sumex
    im[t,-1] = sumim
    exi[t,-1] = sumex_i
    exf[t,-1] = sumex_f
    imi[t,-1] = sumim_i
    imf[t,-1] = sumim_f
    gtrd[t,-1] = sumgtrd
    itrd[t,-1] = sumitrd
    ftrd[t,-1] = sumftrd
    ishare[t,-1] = sumitrd
    tbe[t,-1] = 0

data={'go':go,
      'gdp':gdp,
      'ex':ex,
      'im':im,
      'exi':exi,
      'exf':exf,
      'imi':imi,
      'imf':imf,
      'tbe':tbe,
      'gtrd':gtrd,
      'itrd':itrd,
      'ftrd':ftrd,
      'ishare':ishare}

with open(outpath+'proddata_regions.txt','wb') as f:
    s='region,year'
    for key in data.keys():
        s = s+','+key
    f.write(s+'\n')

    for i in range(0,nr+1):
        for t in range(len(years)):
            s=regions_[i]+','+str(years[t])
            for val in data.values():
                s=s+','+str(val[t,i])
            f.write(s+'\n')

for key,val in data.iteritems():
    with open(outpath+key+'-regions.txt','wb') as f:
        s='Year'
        for r in regions_:
            s = s+','+r
        f.write(s+'\n')

        for t in range(len(years)):
            s=str(years[t])
            for i in range(nr+1):
                s = s+','+str(val[t,i])
            f.write(s+'\n')
