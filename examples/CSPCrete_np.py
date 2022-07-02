# This is the new version of the CSPCret.py
#%%
import pathlib
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import make_interp_spline, BSpline

from CombiCSP.misc import  heatmap2d
from CombiCSP import SolarTroughCalcs, SolarTowerCalcs, HOYS_DEFAULT
import CombiCSP.SolarGeometry as sgh
from CombiCSP.storage import Tr

# from Demand_supply import *
#%% Load data and constants
hoy = HOYS_DEFAULT
# read data from local file
"tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
df_pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
Ib = df_pvgis.loc[:,'Gb(n)']

#%% Tower related dimensions
# Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
# alt = 200*10e-3 #Height above sea level [m] # TODO this is probably 200*1e-3
# Ht = 0.1 #np.arange(0.1,0.4,0.1) # Tower height [km]
# A_helio = 225000 # SolarII 82,750 m² for 10MW https://en.wikipedia.org/wiki/The_Solar_Project

stc =  SolarTowerCalcs(alt = 200*10e-3 , Ht = 0.1, Ar = 99.3 , A_helio = 225000)
oTow = stc.perform_calc(Ib)
#%%
plt.plot(hoy, oTow.data, label='1')
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.title('Tower')
plt.legend()


#%%
# Trough dimensions
# foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
# Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
# L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine
# N = 1800 # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine
# Wr_ini=0.07 # tube outer diameter [m]
# Wc_ini=5.76 # collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM


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
plt.plot(hoy, ons.data)#,xlim(100,600)
plt.plot(hoy, oew.data)#,xlim(100,600)
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.legend(('EW','NS'))
plt.title('Trough')
#xlim(0,87.60), ylim(0,80)
#np.savetxt('datah.txt',datah.T ,delimiter=',') #save transposed data
#np.savetxt('dataxyz.txt',datah.T ,delimiter=',') #save tro_xyz data
plt.show()



plt.title('Tower')
heatmap2d(oTow.data4surf().T)
plt.title('Trough N-S')
heatmap2d(ons.data4surf().T)
plt.title('Trough E-W')
heatmap2d(oew.data4surf().T)


# # %%


# CombiCSP.misc.heatmap_sns(tro_xyz.T, title='Trough E-W')
# # %%
# %%
