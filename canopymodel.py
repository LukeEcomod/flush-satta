# -*- coding: utf-8 -*-
"""
Created on Fri Apr 28 16:03:39 2017

@author: slauniai


*****************
MAIN CODE FOR RUNNING CANOPYGRID FOR SATTASUO. Called by Sattasuo_main

NEEDS FOLLOWING PACKAGES: numpy, pandas, os, matplotlib, configparser

"""


import numpy as np
import pandas as pd
import configparser
import timeit
# from Spathy_utils import initialize_netCDF

import CanopyGrid_satta as cgrid # model core, modified 
    
eps = np.finfo(float).eps  # machine epsilon


def run_model(setupfile):
    
    """ CanopyGrid run for Sattasuo """

    #setupfile = 'c:\\projects\\FLUSH-peatland\\Sattasuo\\model\\sattasuo.ini'
    
    # following creates 1 x 1 grid with constant parameters read from file 
    cmask = np.ones(1)
    ba = None
    lai = None
    hc = None    
    laiu = None
    # for gridded simulations, give here:
    # cmask = sattasuo comp. grid with ones
    # lai = local leaf mass / m2 x SLA
    # ba = local basal area
    # hc = local canopy height but since hc is not used yet give just constant value in grid with shape(cmask)
    
    start_time = timeit.default_timer()
    
    """read ini-file and create model instance 'run1' """
    pgen, pcpy = read_setup(setupfile)
    dt = pgen['dt']
    
    run1 = cgrid.CanopyGrid(pcpy, lai, laiu, ba, hc,  cmask)  

    """ read FMI 10x10grid data for sattasuo"""
    FORC =  read_FMI_weatherdata(pgen['catchment_id'], pgen['start_date'],pgen['end_date'], pgen['forcing_file'])
    FORC['Prec']=FORC['Prec']/dt #mms-1
    Nsteps = len(FORC['Prec'])  # number of timesteps


#    """ create netCDF output file """
#    ncf,outf=initialize_netCDF(spa.id,spa.GisData,spa.FORC, fpath=spa.pgen['output_folder']); #netCDF-file handle
        
    """ initialize dictionary for outputs """

    res = {'time': FORC.index, 'SWE': np.ones(Nsteps)*np.NaN, 'PotInf': np.ones(Nsteps)*np.NaN, 
           'Trfall': np.ones(Nsteps)*np.NaN, 'Interc_o': np.ones(Nsteps) * np.NaN, 
           'Interc_u': np.ones(Nsteps) * np.NaN,'ET': np.ones(Nsteps)*np.NaN,
           'Tr_o': np.ones(Nsteps), 'Tr_u': np.ones(Nsteps)*np.NaN, 'Evap_o': np.ones(Nsteps),
           'Evap_u': np.ones(Nsteps)*np.NaN, 'Mbe': np.ones(Nsteps)*np.NaN,
           'W': np.ones(Nsteps)*np.NaN, 'Wu': np.ones(Nsteps)*np.NaN
           }
     
    
    print '******* Running Sattasuo simulations ********'
    
    for k in range(Nsteps):
        #print 'k='+str(k)

        #forcing
        doy = FORC['doy'].iloc[k]
        ta = FORC['T'].iloc[k]
        vpd = FORC['VPD'].iloc[k]
        rg = FORC['Rg'].iloc[k]
        par = FORC['Par'].iloc[k]
        prec = FORC['Prec'].iloc[k]
        Rew0 = 1.0  # relatively extractable water; set to 1 omits soil restrictions to transpi        

        #run CanopyGrid
        potinf, trfall, interc_o, interc_u, evap_o, evap_u, et, tr_o, tr_u,  mbe = \
            run1.run_CanopyGrid(doy, dt, ta, prec, rg, par, vpd, Rew=Rew0, P=101300.0)

        #save results
        res['SWE'][k] = run1.SWE
        res['W'][k] = run1.W
        res['Wu'][k] = run1.Wu
        res['PotInf'][k] = potinf
        res['Trfall'][k] = trfall            
        res['Interc_o'][k] = interc_o
        res['Interc_u'][k] = interc_u
        res['ET'][k] = et
        #res['Efloor'][k] = ef
        res['Tr_o'][k] = tr_o
        res['Tr_u'][k] = tr_u
        res['Evap_o'][k] = evap_o
        res['Evap_u'][k] = evap_u
        res['Mbe'][k] = mbe
                       
    print('Loops total [s]: ',(timeit.default_timer() - start_time)  )
    print '********* done *********'
    
    return run1, res, FORC


""" ******** Utility functions **************** """

    
def read_setup(inifile):
    """read_setup(inifile): reads .ini parameter file into pp dict"""
    
    cfg = configparser.ConfigParser()
    cfg.read(inifile)

    pp = {}
    for s in cfg.sections():
        section = s.encode('ascii', 'ignore')
        pp[section] = {}
        for k, v in cfg.items(section):
            key = k.encode('ascii', 'ignore')
            val = v.encode('ascii', 'ignore')
            if section == 'General':
                pp[section][key] = val
            else:
                pp[section][key] = float(val)

    pp['General']['dt'] = float(pp['General']['dt'])

    pgen = pp['General']
    pcpy = pp['CanopyGrid']

    return pgen, pcpy


""" ************************ Forcing data, sitefile ************************** """


def read_FMI_weatherdata(ID, start_date, end_date, sourcefile='Satta_FMI_10x10km.csv'):
    """ 
    reads FMI interpolated daily weather data from file 
    IN: 
        ID - sve catchment ID. set ID=0 if all data wanted
        start_date - 'yyyy-mm-dd'
        end_date - 'yyyy-mm-dd'
    OUT:
        fmi - pd.dataframe with datetimeindex
            fmi columns:['ID','Kunta','aika','lon','lat','T','Tmax','Tmin','Prec','Rg','h2o','dds','Prec_a','Par','RH','esa','VPD','doy']
            units: T, Tmin, Tmax, dds[degC], VPD, h2o,esa[kPa], Prec, Prec_a[mm], Rg,Par[Wm-2],lon,lat[deg]
    """
    
    #OmaTunniste;OmaItä;OmaPohjoinen;Kunta;siteid;vuosi;kk;paiva;longitude;latitude;t_mean;t_max;t_min;
    #rainfall;radiation;hpa;lamposumma_v;rainfall_v;lamposumma;lamposumma_cum
    #-site number
    #-date (yyyy mm dd)
    #-latitude (in KKJ coordinates, metres)
    #-longitude (in KKJ coordinates, metres)
    #-T_mean (degrees celcius)
    #-T_max (degrees celcius)
    #-T_min (degrees celcius)
    #-rainfall (mm)
    #-global radiation (per day in kJ/m2)
    #-H2O partial pressure (hPa)

    ID=int(ID)
        
    #import forcing data
    fmi = pd.read_csv(sourcefile, sep=';', header='infer', usecols=['OmaTunniste','Kunta','aika','longitude','latitude','t_mean','t_max','t_min',\
    'rainfall','radiation','hpa','lamposumma_v','rainfall_v'],parse_dates=['aika'])
    
    time = pd.to_datetime(fmi['aika'],format='%Y%m%d')
    
    fmi.index = time
    fmi = fmi.rename(columns={'OmaTunniste': 'ID', 'longitude':'lon','latitude':'lat','t_mean':'T','t_max':'Tmax','t_min':'Tmin','rainfall':'Prec',\
        'radiation':'Rg', 'hpa':'h2o','lamposumma_v':'dds', 'rainfall_v': 'Prec_a'})
    
    fmi['h2o']=1e-1*fmi['h2o'] #hPa-->kPa
    fmi['Rg']=1e3/86400.0*fmi['Rg'] #kJ/m2/d-1 to Wm-2 
    fmi['Par']=0.5*fmi['Rg']

    #saturated vapor pressure    
    esa=0.6112*np.exp((17.67*fmi['T'])/ (fmi['T'] +273.16 -29.66))  #kPa
    vpd=esa - fmi['h2o']; #kPa   
    vpd[vpd<0]=0.0
    rh=100.0*fmi['h2o']/esa;
    rh[rh<0]=0.0; rh[rh>100]=100.0
                
    fmi['RH']=rh;
    fmi['esa']=esa;
    fmi['VPD']=vpd

    fmi['doy']=fmi.index.dayofyear
    #fmi=fmi.drop(['aika'])
    #replace nan's in prec with 0.0
    fmi['Prec'][np.isnan(fmi['Prec'])]=0.0
    #del dat, fields, n, k, time
    
    #get desired period
    fmi=fmi[(fmi.index >= start_date) & (fmi.index <= end_date)]
    if ID >0:
        fmi=fmi[fmi['ID']==ID]
    return fmi
