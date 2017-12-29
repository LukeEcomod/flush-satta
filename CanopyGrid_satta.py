# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 11:01:50 2017

@author: slauniai

******************************************************************************
    SpatHy aboveground domain adapded for Sattasuo -case
    
******************************************************************************

"""
import numpy as np
import configparser
eps = np.finfo(float).eps


class CanopyGrid():
    def __init__(self, cpara, lai=None, laiu = None, ba=None, hc=None,  cmask=None):
        """
        initializes CanopyGrid -object
        IN:
            cpara - parameter dict:
                {'spatial': no, 'wmaxsnow': 1.0, 'lon': 24.21, 'cf': 0.7,
                'wmax': 0.5, 'albedo': 0.1, 'swe': 40.0,pmax': 20.0,
                'q50': 30.0, 'lai': 4.0, 'kfreeze': 5.79e-06, 'elev': 100.0,
                'lat': 60.38, 'hc': 22.0, 'dt': 86400.0, 'kmelt': 2.8934e-05,
                'ka': 0.6, 'f': 0.4, 'emi': 0.98, 'kp': 0.6, 'gsref': 0.0018, 'gsrefu': 0.004,
                'w': 0.0, 'clump': 0.7}
            lai - leaf area index grid; if spatial: yes
            ba - basal area grid
            hc - canopy height grid
            cmask - catchment mask (when lai=None)
        OUT:
            self - object
        self.
            LAI - leaf area index [m2m-2]
            hc - canopy heigth [m]
            cf - closure [-]
            zo - roughness length [m]
            d - displacement height [m]
            clump - clumping index [-]
            physpara - physiologic parameters dict
            phenopara - phenology param. dict
            cpara - copy of inputs, for testing snow model
            Wmax - maximum water storage capacity [mm = kgm-2(ground)]
            WmaxSnow -  maximum snow storage capacity [mm]
            Kmelt - degree day snowmelt factor [mm d-1 K-1]
            Kfreeze - degree day freezing factor [mm d-1 K-1]
            R - snow water retention coeff [-]
            
            State variables:

            X - 'phenological state', i.e. delayed temperature [degC]
            W - canopy water or snow storage [mm]
            SWE - snow water equivalent [mm]
            SWEi - snow water as ice [mm]
            SWEl - snow water as liquid [mm]
            DDsum - degree-day sum [degC]
        """
        epsi = 0.01

        if lai is not None:
            # print '**** CanopyGrid - stands from lai data *******'
            self.LAI = lai + epsi  # m2m-2
            self.hc = hc + epsi
            self.cf = 0.1939 * ba / (0.1939 * ba + 1.69) + epsi
            # canopy closure [-] as function of basal area ba m2ha-1;
            # fitted to Korhonen et al. 2007Silva Fig.2
            
            self.LAIu = laiu + epsi
        else:  # spatially constant properties, given in cpara
            # print '**** CanopyGrid - stands with constant values *******'
            self.LAI = cpara['lai'] * cmask
            self.hc = cpara['hc'] * cmask
            self.cf = cpara['cf'] * cmask
            self.LAIu = cpara['laiu'] * cmask
            
        gridshape = np.shape(self.LAI)
        # print gridshape    
        self.zo = 1.0 / 30.0 * self.hc  # m, roughness length
        self.d = 0.67 * self.hc  # m, displacement height
        
        self.Clump = cpara['clump']  # clumping index
        
        # physiology: transpi + floor evap
        self.physpara = {'q50': cpara['q50'],
                         'gsref': cpara['gsref'], 'gsrefu': cpara['gsrefu'],
                         'ka': cpara['ka'], 'kp': cpara['kp'], 'f': cpara['f'],
                         'rw': cpara['rw'],'rwmin': cpara['rwmin'], 'ga': cpara['ga']}
       
        self.phenopara = {'Smax': cpara['smax'], 'Xo': cpara['xo'],
                          'tau': cpara['tau'], 'fmin': cpara['fmin']}
        
        self.cpara = cpara
        self.Wmax = cpara['wmax'] * self.LAI
        self.Wmaxsnow = cpara['wmaxsnow'] * self.LAI
        self.Wmax_under = cpara['wmax'] * self.LAIu + cpara['wmaxmoss']

        self.Kmelt = cpara['kmelt']
        self.Kfreeze = cpara['kfreeze']
        self.R = 0.05  # max fraction of liquid water in snow

        # --- state variables
        self.X = 0.0
        self.W = np.minimum(cpara['w'], self.Wmax)
        self.Wu = 0.0  # understory storage mm
        self.SWE = cpara['swe']*np.zeros(gridshape)
        self.SWEi = self.SWE
        self.SWEl = np.zeros(gridshape)
        self.DDsum = 0.0

            
    def run_CanopyGrid(self, doy, dt, Ta, Prec, Rg, Par, VPD, Rew=1.0,
                       P=101300.0):
        """
        Runs CanopyGrid instance for one timestep
        IN:
            doy - day of year
            dt - timestep [s]
            Ta - air temperature  [degC], scalar or (n x m) -matrix
            prec - precipitatation rate [mm/s]
            Rg - global radiation [Wm-2], scalar or matrix
            Par - photos. act. radiation [Wm-2], scalar or matrix
            VPD - vapor pressure deficit [kPa], scalar or matrix
            Rew - relative extractable water [-], scalar or matrix
            P - pressure [Pa], scalar or matrix
        OUT:
            updated CanopyGrid instance state variables
            flux grids PotInf, Trfall, Interc, Evap, ET, MBE [mm]
        """
        Cair = P / (8.31 * (Ta + 273.15))  # molar concentration of air, mol m-3
        Mw = 18e-3  #  H2O molar mass [kg mol-2]
        C = Cair * Mw  # H2O flux from ms-1 to kg m-2 s-1 = mm s-1
        
        Rn = 0.7 * Rg #net radiation: replace with function of LAI
#        Rn = np.maximum(2.57 * self.LAI / (2.57 * self.LAI + 0.57) - 0.2,
#                        0.55) * Rg  # Launiainen et al. 2016 GCB Fig 2a fit

        Ga = self.physpara['ga']  # aerodynamic conductance, assume fixed [ms-1]  
        # f = self.physpara['f'] * np.ones(np.shape(self.SWE))
        # f[self.SWE > 0] = eps  # in case of snow, neglect forest floor evap

        """ --- phenology: update self.DDsum & self.X ---"""
        fPheno = self._photoacclim(Ta)
        # print doy,self.X, fPheno

        """---interception and snowpack storage---"""

        # amount of net radiation above understory & at forest floor
        Rnu = Rn * np.exp(- self.physpara['kp']*self.Clump * self.LAI)
        #Rnf = Rnu * np.exp(- self.physpara['kp']*self.LAIu)

        # equilibrium evaporation rates from overstory & understory        
        erate_o = 1.26 * eq_evap(Rn - Rnu, Ta, units='mm') * dt  # mm
        erate_u = 1.26 * eq_evap(Rnu, Ta, units='mm') * dt  # mm
        
        erate_o = np.maximum(0.0, erate_o)
        erate_u = np.maximum(0.0, erate_u)
        erate_u[self.SWE > 0] = 0.0  # in case of snow
        # if Ta < 0 : erate = np.minimum(0.0, erate)
        
        # Call function to return:
        #   PotInf - potential infiltation below bottom layer [mm]
        #   Trfall - throughfall above understory [mm]
        #   Evap_ - evaporation from overstory & understory [mm]
        #   Interc_ - interception, over- & understory [mm]
        #   MBE - mass balance error [mm]
        
        PotInf, Trfall, Evap_o, Evap_u, Interc_o, Interc_u, MBE = \
            self.canopy_water_snow(dt, Ta, Prec, erate_o, erate_u)

        del erate_o, erate_u
        
        """
        Conductances [m s-1] & transpiration of canopy layers [mm]
        """
        Gc_over, Gc_under = self.canopy_conductance(VPD, Par, Ga, CO2=380.0, Rew=Rew, fPheno=fPheno)
                
        if Ta > 0:
             tr_o = C *Gc_over * 1e3* VPD / P  * dt  # mm / dt
             tr_u = C * Gc_under * 1e3* VPD / P * dt
        else:
            tr_o = np.zeros(np.shape(self.LAI))
            tr_u = np.zeros(np.shape(self.LAI))

#       # Evaporation from forest floor: this is now removed since moss layer is handled as
#       # interception storage. Assume that living mosses block evaporation from deeper soil; this
#       # is consistent with FLUSH -simulations that 'kill' soil evaporation
#        if self.SWE == 0:
#            Ef = dt *eq_evap(Rnf, Ta, units='mm')
#        else:
#            Ef = np.zeros(np.shape(self.LAI))

        # evapotranspiration [mm]    
        ET = tr_o + tr_u + Evap_o + Evap_u
        
        return PotInf, Trfall, Interc_o, Interc_u, Evap_o, Evap_u, ET, tr_o, tr_u,  MBE


    
    def _photoacclim(self, T):
        """
        computes new stage of temperature acclimation and phenology modifier.
        Peltoniemi et al. 2015 Bor.Env.Res.
        IN: object, T = daily mean air temperature
        OUT: fPheno - phenology modifier [0...1], updates object state                    
        """
          
        self.X = self.X + 1.0 / self.phenopara['tau'] * (T - self.X)  # degC
        S = np.maximum(self.X - self.phenopara['Xo'], 0.0)
        fPheno = np.maximum(self.phenopara['fmin'],
                            np.minimum(S / self.phenopara['Smax'], 1.0))
        return fPheno


    def canopy_conductance(self, D, Qp, Ga, CO2=380.0, Rew=1.0, fPheno=1.0):
        """                
        Computes integrated canopy conductance Gc for vegetation layers 'overstory' and
        'understory'
        IN:
           D - vpd in kPa
           Qp - PAR in Wm-2
           Ga - aerodynamic conductance (ms-1)
           CO2 - co2 mixing ratio (ppm)
           Rew - relative extractable water [-]
           fPheno - phenology modifier [-]
        OUT:
           Gc - overstory canopy conductance (m/s)
           Gcu - understory canopy conductance (m/s)
        NOTE:  Gc is Canopy condutance (integrated stomatal conductance)
        
        SOURCES:
        Launiainen et al. (2016). Do the energy fluxes and surface conductance
        of boreal coniferous forests in Europe scale with leaf area?
        Global Change Biol.
        Modified from: Leuning et al. 2008. A Simple surface conductance model
        to estimate regional evaporation using MODIS leaf area index and the
        Penman-Montheith equation. Water. Resources. Res., 44, W10419
        Original idea Kelliher et al. (1995). Maximum conductances for
        evaporation from global vegetation types. Agric. For. Met 85, 135-147
        
        Samuli Launiainen, Luke 6-11/2015. Converted to Python 21.9.2016
        
        """

        # ---- model parameters:
        gsref = self.physpara['gsref']  # maximum (leaf-scale) gs (m/s) 
        gsrefu = self.physpara['gsrefu']  # maximum understory (leaf-scale) gs (m/s) 
        kp = self.physpara['kp']  # (-) attenuation coefficient for PAR
        q50 = self.physpara['q50']  # Wm-2, half-sat. of leaf light response
        rw = self.physpara['rw']  # rew parameter
        rwmin = self.physpara['rwmin']  # rew parameter
      
        """--- canopy conductance Gc (integrated stomatal conductance)----- """
        # light responses
        # fQ: Saugier & Katerji, 1991 Agric. For. Met., eq. 4. Leaf light response = Qp / (Qp + q50)
        # overstory        
        fQo = 1./ kp * np.log((kp*Qp + q50) / (kp*Qp*np.exp(-kp * self.LAI) + q50 + eps) )
        # understory        
        Qpu = Qp * np.exp(-kp * self.Clump * self.LAI)  # PAR at understory top
        fQu = 1./ kp * np.log((kp*Qpu + q50) / (kp*Qpu*np.exp(-kp * self.LAIu) + q50 + eps) )
        # VPD -response
        fD = 1.0 / (np.sqrt(D) + eps)

        # soil moisture response: Lagergren & Lindroth, xxxx"""
        fRew = np.minimum(1.0, np.maximum(Rew / rw, rwmin))
        # fRew = 1.0

        # CO2 -response of canopy conductance, derived from APES-simulations
        # (Launiainen et al. 2016, Global Change Biology). relative to 380 ppm
        fCO2 = 1.0 - 0.387 * np.log(CO2 / 380.0)
        # fCO2 = 1.0
        
        # conductances for over- and understory [ms-1]
        Gc = gsref * fQo * fD * fRew * fPheno * fCO2
        Gcu = gsrefu * fQu * fD * fRew * fPheno * fCO2

        return Gc, Gcu

    def canopy_water_snow(self, dt, T, Prec, erate_o, erate_u):
        """
        Calculates canopy and field layer water interception and SWE during timestep dt
        INPUT: self - object
               dt - timestep [s]
               T - air temperature (degC)
               Prec - precipitation rate during (mm d-1)
               erate_o - evaporative demand overstory (mm d-1)
               erate_u - evaporative demand field layer (mm d-1)
        OUTPUT:
            self - updated state W, Wf, SWE, SWEi, SWEl
            PotInf - potential infiltration to soil profile (mm)
            Trfall - throughfall above understory (mm)
            Evap_o - evaporation from canopy store (mm)
            Evap_u - -"- from field+bottom layer storage (mm)
            Interc_o - interception at canopy (mm)
            Interc_u - interception at field+bottom layer storage (mm)
            MBE - mass balance error (mm)
        """

        Kmelt = self.Kmelt - 1.64 * self.cf / dt  # Kuusisto Esko, 'Lumi Suomessa -juttu'

        Kfreeze = self.Kfreeze
        Tmin = 0.0  # 'C, below all is snow
        Tmax = 1.0  # 'C, above all is water
        Tmelt = 0.0  # 'C, T when melting starts

        gridshape = np.shape(self.LAI)  # rows,cols        
        Prec = Prec * dt  # mm
        if Prec > 0:
            erate_o = 0.0  # neglect evaporation during precipitation
            erate_u = 0.0

        # ---Estimate state of precipitation [as water (fW) or as snow(fS)]    
        if T >= Tmax:
            fW = 1; fS = 0
            # Wmax=self.Wmax
        elif T < Tmin:
            fW = 0; fS = 1
            # Wmax=self.Wmaxsnow
        else:
            fW = (T - Tmin) / (Tmax - Tmin)
            fS = 1 - fW
            # Wmax = fW * self.Wmax + fS * self.Wmaxsnow
        
        # ---- Local fluxes (mm)
        Unload = np.zeros(gridshape)  # snow unloading
        Interc = np.zeros(gridshape)  # interception
        Interc_u = np.zeros(gridshape)  # interception
        Melt = np.zeros(gridshape)   # melting
        Freeze = np.zeros(gridshape)  # freezing
        Evap_o = np.zeros(gridshape)
        Evap_u = np.zeros(gridshape)        

    
        """ --- initial conditions for calculating mass balance error --"""
        Wo = self.W  # canopy storage mm
        SWEo = self.SWE  # Snow water equivalent mm
        Wuo = self.Wu # field + bottom layer storage mm
        
        """ --------- Overstory canopy water storage ----------------------------"""
        
        # snow unloading from canopy if T>Tmax # and mean T >0 at previous day
        if T > Tmax: 
            Unload = np.maximum(self.W - self.Wmax, 0.0)
            self.W = self.W - Unload

        # Interception of rain or snow: asymptotic approach of saturation.
        # Hedstrom &Pomeroy 1998. Hydrol. Proc 12, 1611-1625;
        # Koivusalo & Kokkonen 2002 J.Hydrol. 262, 145-164.

        if T < Tmin:  # all Prec is snow
            Interc = (self.Wmaxsnow - self.W) \
                    * (1.0 - np.exp(-self.cf / self.Wmaxsnow * Prec))
        else:
            Interc = np.maximum(0.0, (self.Wmax - self.W)
                    * (1.0 - np.exp(-self.cf / self.Wmax * Prec)))
                    
        self.W = self.W + Interc  # new canopy storage, mm
        
        Trfall = Prec + Unload - Interc  # to virtual snowpack, mm
        
    
        # evaporate from storage
        Evap_o = np.minimum(erate_o, self.W)  # mm
        self.W = self.W - Evap_o  # new storage after evaporation

        """ -- Snow processes (in case no snow, all Trfall routed as Trfall_u to field layer) ---
        """
        if T >= Tmelt:
            Melt = np.minimum(self.SWEi, Kmelt * dt * (T - Tmelt))  # mm
        elif T < Tmelt:
            Freeze = np.minimum(self.SWEl, Kfreeze * dt * (Tmelt - T))  # mm

        Sice = np.maximum(0.0, self.SWEi + fS * Trfall + Freeze - Melt)
        Sliq = np.maximum(0.0, self.SWEl + fW * Trfall - Freeze + Melt)
     
        Trfall_u = np.maximum(0.0, Sliq - Sice * self.R)  # mm, to field layer
        Sliq = np.maximum(0.0, Sliq - Trfall_u)  # mm, liquid water in snow
        
        # update Snowpack state variables
        self.SWEl = Sliq
        self.SWEi = Sice
        self.SWE = self.SWEl + self.SWEi

        """ check interception at field layer (understory + bottom layer) """
        Interc_u = np.maximum(0.0, (self.Wmax_under - self.Wu)
                    * (1.0 - np.exp(-1.0 / self.Wmax_under * Trfall_u)))
        self.Wu = self.Wu + Interc_u  # new storage, mm
        PotInf = Trfall_u - Interc_u  # potential infiltration to soil [mm]
        
        # evaporate from storage
        Evap_u = np.minimum(erate_u, self.Wu)  # mm
        if self.SWE > 0:
            Evap_u = 0.0
        self.Wu = self.Wu - Evap_u  # new storage after evaporation

        """ mass balance error """
        MBE = (self.W + self.Wu + self.SWE) - (Wo + Wuo + SWEo) - (Prec - Evap_o - Evap_u - PotInf)

        if np.max(MBE) > 0.1:
            print('canopy_water MBE [mm] !', MBE)

        return PotInf, Trfall, Evap_o, Evap_u, Interc, Interc_u, MBE


def read_ini(inifile):
    """read_ini(inifile): reads canopygrid.ini parameter file into pp dict"""
    
    cfg = configparser.ConfigParser()
    cfg.read(inifile)

    pp = {}
    for s in cfg.sections():
        section = s.encode('ascii', 'ignore')
        pp[section] = {}
        for k, v in cfg.items(section):
            key = k.encode('ascii', 'ignore')
            val = v.encode('ascii', 'ignore')
            if section == 'General':  # 'general' section
                pp[section][key] = val
            else:
                pp[section][key] = float(val)
                
    pp['General']['dt'] = float(pp['General']['dt'])

    pgen = pp['General']
    cpara = pp['CanopyGrid']
    return pgen, cpara

# @staticmethod
def degreeDays(dd0, T, Tbase, doy):
    """
    Calculates degree-day sum from the current mean Tair.
    INPUT:
        dd0 - previous degree-day sum (degC)
        T - daily mean temperature (degC)
        Tbase - base temperature at which accumulation starts (degC)
        doy - day of year 1...366 (integer)
    OUTPUT:
        x - degree-day sum (degC)
   """
   
    if doy == 1:  # reset in the beginning of the year
        dd0 = 0.0  
    return dd0 + max(0, T - Tbase)


# @staticmethod
def eq_evap(AE, T, P=101300.0, units='W'):
    """
    Calculates the equilibrium evaporation according to 
    McNaughton & Spriggs, 1986.
    INPUT:
        AE - Available energy (Wm-2)
        T - air temperature (degC)
        P - pressure (Pa)
        units - W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
    OUTPUT:
        equilibrium evaporation rate (Wm-2)
        
        NT=273.15; %0 degC in K
    """
    Mw = 18e-3  # kg mol-1
    # latent heat of vaporization of water [J/kg]
    L = 1e3 * (2500.8 - 2.36 * T + 1.6e-3 * T ** 2 - 6e-5 * T ** 3)
    # "" sublimation [J/kg]    
    if T < 0:
        L = 1e3 * (2834.1 - 0.29 * T - 0.004 * T ** 2)

    _, s, g = e_sat(T, P)
    
    x = np.divide((AE * s), (s + g))  # Wm-2 = Js-1m-2
    if units == 'mm':
        x = x / L  # kg m-2 s-1 = mm s-1
    elif units == 'mol':
        x = x / L / Mw  # mol m-2 s-1
    x = np.maximum(x, 0.0)
    return x
                                         

# @staticmethod    
def e_sat(T, P=101300):
    """
    Computes saturation vapor pressure (Pa), slope of vapor pressure curve 
    [Pa K-1]  and psychrometric constant [Pa K-1]
    IN:
        T - air temperature (degC)
        P - ambient pressure (Pa)  
    OUT:
        esa - saturation vapor pressure in Pa
        s - slope of saturation vapor pressure curve (Pa K-1)
        g - psychrometric constant (Pa K-1)
    """    
    NT = 273.15
    cp = 1004.67  # J/kg/K
    
    Lambda = 1e3 * (3147.5 - 2.37 * (T + NT))  # lat heat of vapor [J/kg]
    esa = 1e3 * (0.6112 * np.exp((17.67 * T) / (T + 273.16 - 29.66)))  # Pa

    s = 17.502 * 240.97 * esa / ((240.97 + T) ** 2)
    g = P * cp / (0.622 * Lambda)
    return esa, s, g


# @staticmethod
def penman_monteith(AE, D, T, Gs, Ga, P=101300.0, units='W'):
    """    
    Computes latent heat flux LE (Wm-2) i.e evapotranspiration rate ET (mm/s) 
    from Penman-Monteith equation
    INPUT:
       AE - available energy [Wm-2]
       VPD - vapor pressure deficit [Pa]
       T - ambient air temperature [degC]
       Gs - surface conductance [ms-1]
       Ga - aerodynamic conductance [ms-1]
       P - ambient pressure [Pa]
       units - W (Wm-2), mm (mms-1=kg m-2 s-1), mol (mol m-2 s-1)
    OUTPUT:
       x - evaporation rate in 'units'
    """
    # 'constants'
    cp = 1004.67  # J kg-1 K-1
    rho = 1.25  # kg m-3
    Mw = 18e-3  # kg mol-1
    _, s, g = e_sat(T, P)  # slope of sat. vapor pressure, psycrom const
    L = 1e3 * (3147.5 - 2.37 * (T + 273.15))

    x = (s * AE + rho * cp * Ga * D) / (s + g * (1.0 + Ga / Gs))  # Wm-2
    
    if units is 'mm':
        x = x / L  # kgm-2s-1 = mms-1
    if units is 'mol':
        x = x / L / Mw  # mol m-2 s-1
    
    x = np.maximum(x, 0.0)
    return x


# @staticmethod
def aerodynamic_conductance_from_ust(Ust, U, Stanton):
    """    
    computes canopy aerodynamic conductance (ms-1) from frict. velocity
    IN:
       Ustar - friction velocity (ms-1)
       U - mean wind speed at flux measurement heigth (ms-1)
       Stanton - Stanton number (kB-1) for quasi-laminar boundary layer
           resistance. Typically kB=1...12, use 2 for vegetation ecosystems 
           (Verma, 1989, Garratt and Hicks, 1973)
    OUT:
       Ga - aerodynamic conductance [ms-1]
    """   
    kv = 0.4  # von Karman constant
    ra = U / (Ust ** 2 + eps) + Stanton / (kv * (Ust + eps))  # sm-1
    Ga = 1.0 / ra  # ms-1
    return Ga
