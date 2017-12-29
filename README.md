# flush-satta
Computing inputs to FLUSH 3D soil model to study hydrology of drained peatland forest (Stenberg et al. 2018)

The repository contains Python-codes, ini-files and input/output datafiles used to produce input data to FLUSH 3D soil model at Sattasuo drained peatland forest.

Satta_main.py is the main code, calls run_model from canopymodel.py.
run_model.py is driver for canopy interception and evapotranspiration model CanopyGrid. The class is modified version of SpatHy -code.

Inputs: daily weather data (Precip, Tair, Rglob, VPD) & parameters defined in sattasuo.ini
* doy
* Prec - mm d-1
* Rg - global radiation Wm-2
* T, Tmax, Tmin - degC mean, max and minimum daily temperature
* VPD - kPa vapor pressure deficit

Calculates:
* interception of rain & snow in overstory canopy
* throughfall rate to understory layer
* virtual snowpack (snow water equivalent SWE, mm) evolution at forest floor
* interception of rain or meltwater at field/bottomlayer
* potential infiltration to soil
* evaporation from free water storages (canopy & field/bottom layer) is basen on Priestley-Taylor equilibrium evaporation; attenuation of radiation in canopy is exponential
* Transpiration from over- and understory layers (Tr_o, Tr_u, mmd-1) is computed assuming well-coupled conditions:

Tr_i = C x Gi x VPD / P, 

where Gi is layer's integrated conductance (ms-1), VPD vapor pressure deficit (kPa) and P (kPa) ambient pressure. C = Cair * Mwater is conversion factor from ms-1 to kgm-2s-1.  

The Gi is based on Launiainen et al. (2018) HESS, modified from Saugier & Katerji, 1991 and Kelliher et al. 1995.

Gi = gsref * fQ * fVPD * fPheno * fRew

gsref = leaf stomatal conductance at 1kPa and saturated light (ms-1)
fQ = light modifier, fVPD = vpd modifier, fPheno = phenology modifier, fRew = soil water modifier (here fRew = 1 since FLUSH accounts for restriction)

Returns:

* PotInf - mm d-1 potential infiltration rate to soil (i.e. throughfall below living moss layer)
* Tr_o - mm d-1 overstory transpiration rate 
* Tr_u - mm d-1 understory transpiration rate
* Trfall - mm d-1 throughfall rate below overstory (comparison for measurements)
* Interc_o - mm d-1 interception rate overstory
* Interc_u - mm d-1 interception rate understory
* Evap_o - mm d-1 evaporation rate from ovestory storage
* Evap_u - mm d-1 evaporation rate from understory storage
* ET - mm d-1 total ET
* SWE - mm snow water equivalent

