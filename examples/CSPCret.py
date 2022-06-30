#%%
import pathlib
from scipy import integrate
import numpy as np
# from pylab import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import make_interp_spline, BSpline


# from SolarGeometry_hoy import *
import CombiCSP.SolarGeometry as sgh
from CombiCSP.CSP import *
import CombiCSP.misc
# from Demand_supply import *
#%%
hoy = sgh.HOYS_DEFAULT
# read data from local file
"tmy_35.010_26.130_2007_2016.csv"#Atherinolakos
FNAME = pathlib.Path('example_data/tmy_35.015_25.755_2005_2020.csv')
pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
Ib = pvgis.loc[:,'Gb(n)']

# Tower dimensions
Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
alt = 200*10e-3 #Height above sea level [m]
Ht = 0.1 #np.arange(0.1,0.4,0.1) # Tower height [km]
for R in np.arange(1,2,1):
    A_helio = 225000 # SolarII 82,750 m² for 10MW https://en.wikipedia.org/wiki/The_Solar_Project
    Ctow = A_helio / Ar
    tower = solarII(Ib,1,IAM_tow(hoy),A_helio,Ar)
    #inter = make_interp_spline(hoy, tower, k=3)
    #hoy_new = np.linspace(hoy.min(), hoy.max(), 87600)
    #tower_inter = inter(hoy_new)
    #plot(hoy_new, tower_inter, label=R)
    plt.plot(hoy, tower, label=R)
    tow_datah = np.vstack((hoy, tower))
    tow_xyz = np.vstack(tower).reshape((365,24)) # reshape 8760,1 to 365,24
    Ptower = np.amax(tower) # used in CSPecon .round(2)
    Etower = integrate.trapz(tower).round(2) # used in CSPecon
    CF_tow = Etower / (8760 * Ptower)#.round(2)
    tow_data = np.vstack((A_helio,Ctow,Ptower,Etower,CF_tow)) # vertical stack
plt.xlabel('Time (hour of year)')
plt.ylabel('Power (MW)')
plt.title('Tower')
plt.legend()
#xlim(0,87.60), ylim(0,80)
plt.show()

#%%
# Trough dimensions
foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
Wr_ini=0.07 # tube outer diameter [m]
Wr_end=0.08
Wr_step=0.01
Wc_ini=5.76 # collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM
Wc_end=5.77
Wc_step=0.01
Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine
N = 1800 # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine

for Wr in np.arange(Wr_ini, Wr_end, Wr_step):
    for Wc in np.arange(Wc_ini, Wc_end, Wc_step):
        area = Ac(Wc, L, N)
        trough = di_sst(Ib,costhetai_NS(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
        troughew = di_sst(Ib,costhetai_EW(),IAM_tro(hoy),Tr, Wc, Wr, Ws, L, N)
        plt.plot(hoy, trough)#,xlim(100,600)
        plt.plot(hoy, troughew)#,xlim(100,600)
        datah = np.vstack((hoy, trough))
        tro_xyz = np.vstack(trough).reshape((365,24)) # reshape 8760,1 to 365,24
        tro_xyzew = np.vstack(troughew).reshape((365,24)) # reshape 8760,1 to 365,24
        Ptrough = np.amax(trough) # used in CSPecon .round(2)
        Ptroughew = np.amax(troughew) # used in CSPecon .round(2)
        Etrough = integrate.trapz(trough).round(2) # used in CSPecon
        Etroughew = integrate.trapz(troughew).round(2) # used in CSPecon
        CF_tro = Etrough / (8760 * Ptrough)#.round(2)
        CF_troew = Etroughew / (8760 * Ptroughew)#.round(2)
        tro_data = np.vstack((Ac(Wc, L, N),Cg_tro(Wc, Wr, L, N),Ptrough,Etrough,CF_tro)) # vertical stack
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



plt.title('Tower')
CombiCSP.misc.heatmap2d(tow_xyz.T)
plt.title('Trough N-S')
CombiCSP.misc.heatmap2d(tro_xyz.T)
plt.title('Trough E-W')
CombiCSP.misc.heatmap2d(tro_xyzew.T)


# %%


CombiCSP.misc.heatmap_sns(tro_xyz.T, title='Trough E-W')
# %%
