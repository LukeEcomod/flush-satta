Samuli Launiainen 29.12.2017

Sattasuo interception & water budget computations for FLUSH -inputs

Results in files with following format, Flush inputs in first columns
time - yyyy-mm-dd
PotInf - mm d-1 potential infiltration rate to soil (i.e. throughfall below living moss layer)
Tr_o - mm d-1 overstory transpiration rate 
Tr_u - mm d-1 understory transpiration rate
Trfall - mm d-1 throughfall rate below overstory (comparison for measurements)
Interc_o - mm d-1 interception rate overstory
Interc_u - mm d-1 interception rate understory
Evap_o - mm d-1 evaporation rate from ovestory storage
Evap_u - mm d-1 evaporation rate from understory storage
ET - mm d-1 total ET
SWE - mm snow water equivalent
Prec - mm d-1
Rg - global radiation Wm-2
T, Tmax, Tmin - degC mean, max and minimum daily temperature
VPD - kPa vapor pressure deficit

***********************************************
satta_fluxes.dat

Sattasuo current 2006 - 2011 growing seasons (2.6 - 30.10)
LAIo = 2.6 m2m-2
LAIu = 1.2 m2 m-2
bottom layer living mosses 0.16 kgm-2, interception capacity 15 g/g-1
LAIu and moss biomasses estimated using M&M biomass equations (eq. 42 & 43).
Convert field layer biomass kgm-2 to LAI by SLA ~150 gm-2 (Kolari et al. 2006 ForEco).
LAIu seems realistic, see also Sonnentag et al. 2007 (AFM)
NOTE: magniture and response of Trfall / P to thinning seems realistic compared to 'sadesuppilot -data'
***** GS waterbudget partitioning *******
('ET/P', 0.86817143349840953)
('Trfall/P', 0.83500289080460566)
('PotInfil/P', 0.74662210550832764)
('Tro/P', 0.39330056621050707)
('Tru/P', 0.11622727600214863)
***** ET components *******
('Tro/ET', 0.45302177776761365)
('Tru/ET', 0.13387595066771055)
('Evap_o/ET', 0.28579811940434113)
('Evap_u/ET', 0.12730415216033583)
('Transpi/PotInf', 0.68244408845322169)

satta_fluxes_thinned.dat

sattasuo thinned 2006 - 2011 growing seasons (2.6 - 30.10)
LAIo = 1.6 m2m-2
LAIu = 1.2 m2 m-2
bottom layer living mosses 0.16 kgm-2, interception capacity 15 g/g-1

***** GS waterbudget partitioning *******
('ET/P', 0.74882110277728842)
('Trfall/P', 0.88125208275672295)
('PotInfil/P', 0.77637326808211782)
('Tro/P', 0.28271588050504237)
('Tru/P', 0.15442742362022949)
***** ET components *******
('Tro/ET', 0.37754796099693611)
('Tru/ET', 0.20622739269429846)
('Evap_o/ET', 0.23961219496932132)
('Evap_u/ET', 0.17661245133944328)
('Transpi/PotInf', 0.56305816041960211)