# -*- coding: utf-8 -*-
"""
Created on Tue Dec 19 14:04:48 2017

@author: slauniai

Main code for Sattasuo-simulations: 

1) reading data, plotting figs and saving datafiles
2) need to change setupfile, resfile - and figure names when running current vs. thinned simulations
3) all parameters in ini-files are not necessary - using old code modified for this.

Calls canopymodel that is 'driver' for CanopyGrid - the model class.

Samuli 9.5. - 29.12.2017
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from canopymodel import run_model
  
eps = np.finfo(float).eps  # machine epsilon

os.chdir('c:\\projects\\FLUSH-peatland\\Sattasuo\\model')
# read sapflow data
ffile = 'Sapflow_daily_topython.txt'
sapflow = pd.read_csv(ffile, sep=';', parse_dates=[[0,1,2]], header='infer')
sapflow.index = sapflow['yyyy_mm_dd']
#sapflow = sapflow.drop('yyyy_mm_dd')

# Canopygrid model ini file
setupfile = 'sattasuo.ini'
# setupfile = 'c:\\projects\\FLUSH-peatland\\Sattasuo\\model\\sattasuo_thinned.ini'

# results file
resfile = 'satta_fluxes.dat'

"""
run model, using daily 10x10km FMI weatherdata
Can't use Apukka -data since there is no radiation or Precip measured
"""

dt = 3600.0 * 24  # s

# call model 
run1, res, forc = run_model(setupfile)  # returns model instance, result dict, forcing data

res = pd.DataFrame.from_dict(res)
res.index = res['time']
res = res.drop(['time'], axis=1)
res['Prec'] = forc['Prec'] * dt

plt.figure()
plt.plot(res['Mbe']); plt.title('mass balance error')

#%%  draw few figures and save results to csv file

# select mask for snow-free conditions 1st of June - 30th Oct
gs_mask = ((res.index>'2006-06-01') & (res.index<'2006-10-30')) | ((res.index>'2007-06-01') & (res.index<'2007-10-30')) | \
        ((res.index>'2008-06-01') & (res.index<'2008-10-30')) | ((res.index>'2009-06-01') & (res.index<'2009-10-30')) | \
        ((res.index>'2010-06-01') & (res.index<'2010-10-30')) | ((res.index>'2011-06-01') & (res.index<'2011-10-30'))

# partitioning of water budget during gs_mask
pgs = sum(res['Prec'][gs_mask]) # total Precip, mm
tf = sum(res['Trfall'][gs_mask]) # total throughfall below overstory, mm
tinf = sum(res['PotInf'][gs_mask]) # total potential infiltration to soil, mm

tet = sum(res['ET']) # total ET
tevap_o = sum(res['Evap_o'])  # overstory storage evaporation
tevap_u = sum(res['Evap_u'])  # understory storage evaporation
ttr_o = sum(res['Tr_o'])  # overstory transpiration [mm]
ttr_u = sum(res['Tr_u'])  # understory transpiration [mm]

# ********************************************************
print('***** GS waterbudget partitioning *******')
print('ET/P', tet/pgs)
print('Trfall/P', tf / pgs)
print('PotInfil/P', tinf / pgs)

print('Tro/P', ttr_o / pgs)
print('Tru/P', ttr_u / pgs)

print('***** ET components *******')
print('Tro/ET', ttr_o / tet)
print('Tru/ET', ttr_u / tet)
print('Evap_o/ET', tevap_o / tet)
print('Evap_u/ET', tevap_u / tet)

print('Transpi/PotInf', (ttr_o + ttr_u) / tinf)


# plot figures

plt.figure()
plt.plot(res['ET'], 'ro-', label='ET')
plt.plot(res['Tr_o'], 'go-', label='Tr_over')
plt.plot(res['Tr_u'], 'co-', label='Tr_under')
plt.legend()
plt.ylabel('mm/d')
plt.savefig('Satta_ET.png')


plt.figure()
plt.plot(res['Tr_o'], 'go-', label='Tr_over')
plt.plot(res['Tr_u'], 'co-', label='Tr_under')
plt.plot(sapflow['sapf'], 'mo-')
plt.legend()
plt.ylabel('mm/d')
plt.savefig('Satta_Transpi.png')

plt.figure()
plt.subplot(211)
plt.plot(res['PotInf'], 'ro-', label='PotInfil')
plt.plot(res['Tr_o'] + res['Tr_u'], 'bo-', label='Tr_o + Tr_u')
plt.ylabel('mm d-1')
plt.legend()
plt.subplot(212)
plt.plot(np.cumsum(res['PotInf'] - res['Tr_o'] - res['Tr_u']), 'ro')
plt.ylabel('cumulative (Infil - Tr) [mm]')
plt.savefig('Satta_swbal.png')

# subset of forcing data, save to file
aa = forc[['Rg', 'T', 'Tmax', 'Tmin', 'VPD']]

# re-order res columns
res = res[['PotInf', 'Tr_o', 'Tr_u', 'Trfall',  'Interc_o', 'Interc_u', 'Evap_o', 'Evap_u', 'ET', 'SWE', 'Prec']]

# save results: merge dataframes
res = pd.DataFrame.merge(res, aa, left_index=True, right_index=True)
res.to_csv(resfile, sep=';')