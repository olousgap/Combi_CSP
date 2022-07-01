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
from CombiCSP.CSP import *
import CombiCSP.misc
# from Demand_supply import *
#%% Load data and constants
hoy = sgh.HOYS_DEFAULT
# read data from local file
"tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
Ib = pvgis.loc[:,'Gb(n)']

class OutputContainer():
    """this is a container for the data and reshaping them. 

    For simplicity reason this assums that the time index will always be hourly the days of the year

    Returns:
        _type_: _description_
    """    
    hoy = sgh.HOYS_DEFAULT
    def __init__(self, data:pd.Series, A_helio:float, Ctow:float):
        """_summary_

        Args:
            data (np.array): time series
            A_helio (_type_): The heliostats area in m2
            Ctow (_type_): the ratio of 
        """        
        if not data.shape == (8760,):
            raise ValueError("Wrong time series")
        self.data = data
        self.A_helio = A_helio
        self.Ctow = Ctow

    
    def datah(self):
        return np.vstack((self.hoy, self.data))
    
    def tow_xyz(self):
        """stacks the data in days and hours 

        Returns:
            _type_: _description_
        """        
        return np.vstack(self.data).reshape((365,24))
    
    def Ptower(self)->float:
        """returns the maximum power of the Tower

        Returns:
            _type_: _description_
        """        
        return np.amax(self.data) 

    def Etower(self)->float:
        """retunrs the total energy yield of the solar tower. 

        Returns:
           float : the total power per year in [MWh]
        """        
        return integrate.trapz(self.data).round(2)

    def CF_tow(self)->float: 
        return self.Etower() / (8760 * self.Ptower())
    
    def as_df(self)->float: 
        return pd.DataFrame({'Power_MW':self.data}, index= self.hoy)

    def summary_tower_data(self):    
        tow_data = np.vstack((self.A_helio,self.Ctow,self.Ptower(),self.Etower(),self.CF_tow())) # vertical stack
        return tow_data


class SolarTowerCalcs():
    def __init__(self, 
        alt = 200*10e-3 
        , Ht = 0.1
        , Ar = 99.3 
        , A_helio = 225000
        ):
        self.Ar_m2 = Ar# receiver area [m2] pp.44 in Pacheco
        self.alt_m = alt #Height above sea level [m]
        self.Ht_km = Ht # Tower height [km]
        self.A_helio_m2 = A_helio # SolarII 82,750 m² for 10MW https://en.wikipedia.org/wiki/The_Solar_Project
        self.Ctow = self.A_helio_m2 / self.Ar_m2

    def perform_calc(self, Ib, transmittance=1, hoy=sgh.HOYS_DEFAULT):
        self._hourly_results = OutputContainer(data = solarII(Ib=Ib,Trans=transmittance, IAM=cspC.IAM_tow(hoy),
                A_helio=self.A_helio_m2,Ar=self.Ar_m2),
            A_helio=self.A_helio_m2, Ctow=self.Ctow)
        return self._hourly_results

#%%
stc = SolarTowerCalcs(alt = 200*10e-3 , Ht = 0.1, Ar = 99.3 , A_helio = 225000)
oTow = stc.perform_calc(Ib)
#%%

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

plt.plot(hoy, tower, label=R)
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.title('Tower')
plt.legend()
#xlim(0,87.60), ylim(0,80)
plt.show()

#%% Comparison between the results
np.testing.assert_equal(oTow.data.values, desired=tower.values)
np.testing.assert_equal(oTow.tow_xyz(), desired=tow_xyz)
np.testing.assert_equal(oTow.datah(), desired=tow_datah)
assert oTow.Ptower() == Ptower
assert oTow.Etower() == Etower
assert oTow.CF_tow() == CF_tow
assert np.abs(oTow.data.values-tower.values).sum() == 0, 'Discrepancy in the results'
    
#%%

#%%
# Trough dimensions
# foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
# Wr_ini=0.07 # tube outer diameter [m]
# Wr_end=0.08
# Wr_step=0.01
# Wc_ini=5.76 # collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM
# Wc_end=5.77
# Wc_step=0.01
# Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
# L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine
# N = 1800 # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine

# for Wr in np.arange(Wr_ini, Wr_end, Wr_step):
#     for Wc in np.arange(Wc_ini, Wc_end, Wc_step):
#         area = Ac(Wc, L, N)
#         trough = di_sst(Ib,costhetai_NS(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
#         troughew = di_sst(Ib,costhetai_EW(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
#         plt.plot(hoy, trough)#,xlim(100,600)
#         plt.plot(hoy, troughew)#,xlim(100,600)
#         datah = np.vstack((hoy, trough))
#         tro_xyz = np.vstack(trough).reshape((365,24)) # reshape 8760,1 to 365,24
#         tro_xyzew = np.vstack(troughew).reshape((365,24)) # reshape 8760,1 to 365,24
#         Ptrough = np.amax(trough) # used in CSPecon .round(2)
#         Ptroughew = np.amax(troughew) # used in CSPecon .round(2)
#         Etrough = integrate.trapz(trough).round(2) # used in CSPecon
#         Etroughew = integrate.trapz(troughew).round(2) # used in CSPecon
#         CF_tro = Etrough / (8760 * Ptrough)#.round(2)
#         CF_troew = Etroughew / (8760 * Ptroughew)#.round(2)
#         tro_data = np.vstack((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptrough,Etrough,CF_tro)) # vertical stack
#         tro_dataew = np.vstack((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptroughew,Etroughew,CF_troew)) # vertical stack
#         #plot(hoy, plot_features)
# plt.xlabel('Time (hour of year)')
# plt.ylabel('Power (MW)')
# plt.legend(('EW','NS'))
# plt.title('Trough')
# #xlim(0,87.60), ylim(0,80)
# #np.savetxt('datah.txt',datah.T ,delimiter=',') #save transposed data
# #np.savetxt('dataxyz.txt',datah.T ,delimiter=',') #save tro_xyz data
# plt.show()



# plt.title('Tower')
# CombiCSP.misc.heatmap2d(tow_xyz.T)
# plt.title('Trough N-S')
# CombiCSP.misc.heatmap2d(tro_xyz.T)
# plt.title('Trough E-W')
# CombiCSP.misc.heatmap2d(tro_xyzew.T)


# # %%


# CombiCSP.misc.heatmap_sns(tro_xyz.T, title='Trough E-W')
# # %%
