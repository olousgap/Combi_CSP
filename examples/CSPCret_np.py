#%%
import pathlib
from scipy import integrate
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import make_interp_spline, BSpline

# from SolarGeometry_hoy import *
import CombiCSP.SolarGeometry as sgh
import CombiCSP.CSP as cspC
from CombiCSP.CSP import solarII, di_sst,IAM_tow, IAM_tro, costhetai_NS, costhetai_EW, Ac, Cg_tro
from CombiCSP.storage import Tr
from CombiCSP.SolarTrough import SolarTroughCalcs
# import CombiCSP.misc
# TODO remove `from CombiCSP.CSP import *` especially parameter definitions like T
# TODO split CSP to SolarTower, and SolarTrough

# from Demand_supply import *
#%% Load data and constants
hoy = sgh.HOYS_DEFAULT
# read data from local file
"tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
Ib = pvgis.loc[:,'Gb(n)']

#%%
stc = cspC.SolarTowerCalcs(alt = 200*10e-3 , Ht = 0.1, Ar = 99.3 , A_helio = 225000)
oTow = stc.perform_calc(Ib)
#%% This section contains the old code and the comparison

# Tower dimensions 
Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
alt = 200*10e-3 #Height above sea level [m] # TODO this is probably 200*1e-3
Ht = 0.1 #np.arange(0.1,0.4,0.1) # Tower height [km]

A_helio = 225000 # SolarII 82,750 m² for 10MW https://en.wikipedia.org/wiki/The_Solar_Project
Ctow = A_helio / Ar
tower = solarII(Ib,1, IAM_tow(hoy),A_helio,Ar)
#inter = make_interp_spline(hoy, tower, k=3)
#hoy_new = np.linspace(hoy.min(), hoy.max(), 87600)
#tower_inter = inter(hoy_new)
#plot(hoy_new, tower_inter, label=R)
tow_datah = np.vstack((hoy, tower))
tow_xyz = np.vstack(tower).reshape((365,24)) # reshape 8760,1 to 365,24
Ptower = np.amax(tower) # used in CSPecon .round(2)
Etower = integrate.trapz(tower).round(2) # used in CSPecon
CF_tow = Etower / (8760 * Ptower)#.round(2)
tow_data = np.vstack((A_helio,Ctow,Ptower,Etower,CF_tow)) # vertical stack

plt.plot(hoy, tower, label='1')
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.title('Tower')
plt.legend()
#xlim(0,87.60), ylim(0,80)
plt.show()

# Comparison between the results
# if this does not present a problem then the results are ok
np.testing.assert_equal(oTow.data.values, desired=tower.values)
np.testing.assert_equal(oTow.data4surf(), desired=tow_xyz)
np.testing.assert_equal(oTow.hour_power_arr(), desired=tow_datah)
assert oTow.PowerMax_MW == Ptower
assert oTow.Energy_MWh == Etower
assert oTow.CF == CF_tow
assert np.abs(oTow.data.values-tower.values).sum() == 0, 'Discrepancy in the results'
np.testing.assert_equal(oTow.summary_tower_data(), desired=tow_data)
#%%
plt.plot(hoy, oTow.data, label='1')
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.title('Tower')
plt.legend()


#%%
sotr = SolarTroughCalcs(
        foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
        ,N = 1800 # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine
        ,L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine  
        ,Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
        ,Wr = 0.07 # tube outer diameter [m]
        ,Wc = 5.76
        )
oew = sotr.perform_calcs_EW(Ib=Ib, Tr=Tr)
ons = sotr.perform_calcs_NS(Ib=Ib, Tr=Tr)
#%%
# Trough dimensions
foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine
N = 1800 # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine

Wr_ini=0.07 # tube outer diameter [m]
Wr_end=0.08
Wr_step=0.01
Wc_ini=5.76 # collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM
Wc_end=5.77
Wc_step=0.01

for Wr in np.arange(Wr_ini, Wr_end, Wr_step):
    for Wc in np.arange(Wc_ini, Wc_end, Wc_step):
        area = Ac(Wc, L, N)
        trough   = di_sst(Ib,costhetai_NS(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
        troughew = di_sst(Ib,costhetai_EW(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
        plt.plot(hoy, trough)#,xlim(100,600)
        plt.plot(hoy, troughew)#,xlim(100,600)
        #NS
        datah = np.vstack((hoy, trough))
        tro_xyz = np.vstack(trough).reshape((365,24)) # reshape 8760,1 to 365,24
        Ptrough = np.amax(trough) # used in CSPecon .round(2)
        Etrough = integrate.trapz(trough).round(2) # used in CSPecon
        CF_tro = Etrough / (8760 * Ptrough)#.round(2)
        tro_data = np.vstack((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptrough,Etrough,CF_tro)) # vertical stack
        #EW
        tro_xyzew = np.vstack(troughew).reshape((365,24)) # reshape 8760,1 to 365,24
        Ptroughew = np.amax(troughew) # used in CSPecon .round(2)
        Etroughew = integrate.trapz(troughew).round(2) # used in CSPecon
        CF_troew = Etroughew / (8760 * Ptroughew)#.round(2)
        tro_dataew = np.vstack((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptroughew,Etroughew,CF_troew)) # vertical stack
        #plot(hoy, plot_features)

plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.legend(('EW','NS'))
plt.title('Trough')
#xlim(0,87.60), ylim(0,80)
#np.savetxt('datah.txt',datah.T ,delimiter=',') #save transposed data
#np.savetxt('dataxyz.txt',datah.T ,delimiter=',') #save tro_xyz data
plt.show()



# plt.title('Tower')
# CombiCSP.misc.heatmap2d(tow_xyz.T)
# plt.title('Trough N-S')
# CombiCSP.misc.heatmap2d(tro_xyz.T)
# plt.title('Trough E-W')
# CombiCSP.misc.heatmap2d(tro_xyzew.T)


# # %%


# CombiCSP.misc.heatmap_sns(tro_xyz.T, title='Trough E-W')
# # %%

# %%
# Comparison between the results NS
# if this does not present a problem then the results are ok
np.testing.assert_equal(ons.data.values, desired=trough.values)
np.testing.assert_equal(ons.data4surf(), desired=tro_xyz)
np.testing.assert_equal(ons.hour_power_arr(), desired=datah)
assert ons.PowerMax_MW == Ptrough
assert ons.Energy_MWh == Etrough
assert ons.CF == CF_tro
assert np.abs(ons.data.values-trough.values).sum() == 0, 'Discrepancy in the results'
np.testing.assert_equal(ons.summary_tower_data(), desired=tro_data)
#%%

# Comparison between the results EW
# if this does not present a problem then the results are ok
np.testing.assert_equal(oew.data.values, desired=troughew.values)
np.testing.assert_equal(oew.data4surf(), desired=tro_xyzew)
# np.testing.assert_equal(oew.datah(), desired=datah)
assert oew.PowerMax_MW == Ptroughew
assert oew.Energy_MWh == Etroughew
assert oew.CF == CF_troew
assert np.abs(oew.data.values-troughew.values).sum() == 0, 'Discrepancy in the results'
np.testing.assert_equal(oew.summary_tower_data(), desired=tro_dataew)

# %%
