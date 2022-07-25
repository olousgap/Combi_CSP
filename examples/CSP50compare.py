#%%
import pathlib
import pandas as pd
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt

import numpy_financial as npf

from CombiCSP import HOYS_DEFAULT
import CombiCSP.SolarGeometry as sgh
import CombiCSP.misc as cspm
import CombiCSP.economics as cspe

from CombiCSP.solar_tower import solarII,IAM_tow
from CombiCSP.solar_trough import di_sst, IAM_tro, costhetai_NS, costhetai_EW, Ac, Cg_tro

from CombiCSP.storage import Tr

from CSP50_common_econ import *   # this is a development hack to avoid duplication

#import pcm
#%%


hoy = HOYS_DEFAULT
#%% Imported from CSPEcon ==========================================================================================
# from CSPecon import csp_area_costs, power_block_cost, csp_energy_price \
#     ,csp_discount_rate, oil_price, Eoil

# from CSPCret import Ptrough, Etrough



# End of Import from CSPEcon ==========================================================================================
#%%
# read data from local file
"tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
pvgis_data = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') #Atherinolakos
Ib = pvgis_data.loc[:,'Gb(n)']
#Ib = ineichen().dni
capital_csp = 5000000

area_list = []
cash_flow_list = []

# Tower dimensions
Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
alt = 200*10e-3 #Height above sea level [m]
Ht = 0.1 #np.arange(0.1,0.4,0.1) # Tower height [km]
#R = 5

tow_scenaria = []
for A_helio in np.arange(75000,125001,10000): # 100MW np.arange(150000,250001,10000):
    Ctow = A_helio / Ar
    tower = solarII(Ib,1,IAM_tow(hoy),A_helio,Ar)
    #inter = make_interp_spline(hoy, tower, k=3)
    #hoy_new = np.linspace(hoy.min(), hoy.max(), 87600)
    #tower_inter = inter(hoy_new)
    #plot(hoy_new, tower_inter, label=R)
    tow_xyz = np.vstack(tower).reshape((365,24)) # reshape 8760,1 to 365,24
    tow_mon = np.vstack(tower).reshape((365,24)) # reshape 8760,1 to 365,24
    Ptower = np.amax(tower) # used in CSPecon .round(2)
    Etower = integrate.trapz(tower).round(2) # used in CSPecon
    CF_tow = Etower / (8760 * Ptower)#.round(2)
    CF_towh = Etower / (8760 * tower) #<<<<<<<<<<<<<<<<<<<<<<<<<too noisy, test with positive data
    tow_data = np.vstack((A_helio,Ctow,Ptower,Etower,CF_tow)) # vertical stack
    # economics
    capital_csp_tow = A_helio*csp_area_costs + Ptower*power_block_cost
    revenue_csp_tow = cspe.cashflow(Etower,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tow)
    cash_flow_tow = [-capital_csp_tow] + [revenue_csp_tow for i in range(30)]
    dpb_tow = cspe.discounted_payback_period(csp_discount_rate, cash_flow_tow)
    npv_csp_tow = npf.npv(csp_discount_rate, [-capital_csp_tow] + [cspe.cashflow(Etower,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tow) for i in range(30)])
    irr_csp_tow = npf.irr([-capital_csp] + [cspe.cashflow(Etower,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tow) for i in range(30)])
    area_list.append(A_helio)
    cash_flow_list.append(cash_flow_tow)
    tow_scenaria.append((A_helio,Ctow,Ptower,Etower,CF_tow,dpb_tow,npv_csp_tow,irr_csp_tow,cash_flow_tow))
    plt.plot(hoy, tower, label=A_helio)
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)') 
plt.title('Tower')
plt.legend()
#xlim(0,87.60), ylim(0,80)
plt.show()
#%%
# Trough dimensions
foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
Wr = 0.07 # tube outer diameter [m]
Wc = 5.76 # collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM
Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine

trough_scenaria = []
troughew_scenaria = []
for N in np.arange(800,1301,100): # 100MW np.arange(1000,2001,100):
    area = Ac(Wc, L, N)
    trough = di_sst(Ib,costhetai_NS(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
    troughew = di_sst(Ib,costhetai_EW(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
    plt.plot(hoy, trough, label=N)#,xlim(100,600)
    plt.plot(hoy, troughew, label=N)#,xlim(100,600)
    datah = np.vstack((hoy, trough))
    tro_xyz = np.vstack(trough).reshape((365,24)) # reshape 8760,1 to 365,24
    Ptrough = np.amax(trough) # used in CSPecon .round(2)
    Ptroughew = np.amax(troughew) # used in CSPecon .round(2)
    Etrough = integrate.trapz(trough).round(2) # used in CSPecon
    Etroughew = integrate.trapz(troughew).round(2) # used in CSPecon
    CF_tro = Etrough / (8760 * Ptrough)#.round(2)
    CF_troew = Etroughew / (8760 * Ptroughew)#.round(2)
    tro_data = np.vstack((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptrough,Etrough,CF_tro)) # vertical stack
    tro_dataew = np.vstack((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptroughew,Etroughew,CF_troew)) # vertical stack
    # economics
    capital_csp_tro = area*csp_area_costs + Ptrough*power_block_cost
    revenue_csp_tro = cspe.cashflow(Etrough,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tro)
    cash_flow_tro = [-capital_csp_tro] + [revenue_csp_tro for i in range(30)]
    dpb_tro = cspe.discounted_payback_period(csp_discount_rate, cash_flow_tro)
    npv_csp_tro = npf.npv(csp_discount_rate, [-capital_csp_tro] \
        + [cspe.cashflow(Etrough,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tro) for i in range(30)])
    irr_csp_tro = npf.irr([-capital_csp_tro] \
        + [cspe.cashflow(Etrough,csp_energy_price,Eoil,0.4,-oil_price,capital_csp_tro) for i in range(30)])
    area_list.append(area)
    cash_flow_list.append(cash_flow_tro)
    trough_scenaria.append((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptrough,Etrough,CF_tro,dpb_tro,npv_csp_tro,irr_csp_tro,cash_flow_tro))
    troughew_scenaria.append((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptroughew,Etroughew,CF_troew,dpb_tro,npv_csp_tro,irr_csp_tro,cash_flow_tro))
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)'), 
plt.legend()
plt.title('Trough')
#xlim(0,87.60), ylim(0,80)
plt.show()
#%%
DPB = []
for (x,y) in zip(area_list,cash_flow_list):
    DPB.append(cspe.discounted_payback_period(csp_discount_rate, y).round(2))
    #title(z)
    #show()
dpb_table = pd.DataFrame(DPB, index=area_list, columns = ['DPB (years)'])
dpb_table

# T+NS optimum
A_helio_optNS = 125000
N_opt_NS = 800
tower_opt = solarII(Ib,1,IAM_tow(hoy),A_helio_optNS,Ar)
trough_opt = di_sst(Ib,costhetai_NS(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N_opt_NS)
combiNS = tower_opt + trough_opt
combiNS_xyz = np.vstack(combiNS).reshape((365,24)) # reshape 8760,1 to 365,24
plt.title('Tower + Trough N-S')
cspm.heatmap2d(combiNS_xyz.T)

area_combiNS = Ac(Wc, L, N_opt_NS)
PcombiNS = np.amax(combiNS) # used in CSPecon .round(2)
EcombiNS = integrate.trapz(combiNS).round(2) # used in CSPecon

capital_combiNS = (A_helio_optNS+area_combiNS)*csp_area_costs + PcombiNS*power_block_cost
revenue_combiNS = cspe.cashflow(EcombiNS,csp_energy_price,Eoil,0.4,-oil_price,capital_combiNS)
cash_flow_combiNS = [-capital_combiNS] + [revenue_combiNS for i in range(30)]
dpb_combiNS = cspe.discounted_payback_period(csp_discount_rate, cash_flow_combiNS)
npv_combiNS = npf.npv(csp_discount_rate, [-capital_combiNS] \
    + [cspe.cashflow(EcombiNS,csp_energy_price,Eoil,0.4,-oil_price,capital_combiNS) for i in range(30)])
irr_combiNS = npf.irr([-capital_csp] \
    + [cspe.cashflow(EcombiNS,csp_energy_price,Eoil,0.4,-oil_price,capital_combiNS) for i in range(30)])



# T+EW optimum
A_helio_optEW = 125000
N_opt_EW = 800
tower = solarII(Ib,1,IAM_tow(hoy),A_helio_optEW,Ar)
troughew = di_sst(Ib,costhetai_EW(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N_opt_EW)
combiEW = tower + troughew
combiEW_xyz = np.vstack(combiEW).reshape((365,24)) # reshape 8760,1 to 365,24
plt.title('Tower + Trough E-W')
cspm.heatmap2d(combiEW_xyz.T)

area_combiEW = Ac(Wc, L, N_opt_EW)
PcombiEW = np.amax(combiEW) # used in CSPecon .round(2)
EcombiEW = integrate.trapz(combiEW).round(2) # used in CSPecon

capital_combiEW = (A_helio_optEW+area_combiEW)*csp_area_costs + PcombiEW*power_block_cost
revenue_combiEW = cspe.cashflow(EcombiEW,csp_energy_price,Eoil,0.4,-oil_price,capital_combiEW)
cash_flow_combiEW = [-capital_combiEW] + [revenue_combiEW for i in range(30)]
dpb_combiEW = cspe.discounted_payback_period(csp_discount_rate, cash_flow_combiEW)
npv_combiEW = npf.npv(csp_discount_rate, [-capital_combiEW] 
+ [cspe.cashflow(EcombiEW,csp_energy_price,Eoil,0.4,-oil_price,capital_combiEW) for i in range(30)])
irr_combiEW = npf.irr([-capital_csp] 
+ [cspe.cashflow(EcombiEW,csp_energy_price,Eoil,0.4,-oil_price,capital_combiEW) for i in range(30)])

combi_finance = pd.DataFrame((dpb_combiNS,dpb_combiEW,npv_combiNS,npv_combiEW,irr_combiNS,irr_combiEW)).round(2)

# %%
