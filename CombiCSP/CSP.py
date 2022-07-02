#%%
'''Concentrating Solar Power plants                 ALT + SHIFT +0 to unfold levels
                                                    ALT + 0 to fold levels'''
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy import integrate
# from pylab import *

# from SolarGeometry_hoy import *
import CombiCSP.SolarGeometry as sgh
#from pcm import *
#import pcm
from iapws import IAPWS97
#%%
from numpy import log as ln
NaN = np.nan
pi = np.pi
sin = np.sin
cos = np.cos
asin = np.arcsin
acos = np.arccos
tan = np.tan
atan = np.arctan
arctan = np.arctan
atan2 = np.arctan2
rad = np.radians
deg = np.degrees

z= sgh.z
d= sgh.d
W= sgh.W
thetai= sgh.thetai
azim = sgh.azim
ele = sgh.ele


#%%
# generator data
nG = 0.97 # efficiency of generator Mosleh19
SM = 2.5 # Solar Multiplier

def CtoK(c: float|int):
    """convert Celsius to Kelvin

    Args:
        c (float|int): the temperature in degrees Celsius

    Returns:
        float|int: the temperature in Kelvin
    """    
    k = c + 273
    return k
    
def solarII(Ib:pd.Series,Trans:float,IAM:np.array,A_helio:float,Ar:float)->pd.Series:
    """Calculates the power of the solar tower with heliostat
 
    R.K. McGovern, W.J. Smith, Optimal concentration and temperatures of solar thermal power plants,
    Energy Conversion and Management. 60 (2012) 226–232.
    J.E. Pacheco, R.W. Bradshaw, D.B. Dawson, W.D. la Rosa, R. Gilbert, S.H. Goods, P. Jacobs, M.J. Hale, S.A. Jones, G.J. Kolb, M.R. Prairie, H.E. Reilly, S.K. Showalter, L.L. VANT-HULL, 
    Final Test and Evaluation Results from the Solar Two Project, n.d. https://core.ac.uk/reader/193342950 (accessed September 8, 2020).
 


    Args:
        Ib (pd.Series): direct irradiance
        Trans (float): transmissivity 
        IAM (np.array): incidence angle modifier
        A_helio (float): heliostat area in m^2
        Ar (float): receiver area in m^2

    Returns:
        pd.Series: power in MW
    """
    Effopt=100 # heliostat optical effiency [%] 65% pp.24 in Pacheco
    #A_helio = 71140 + 10260 # total heliostat area [m2] pp.22 in Pacheco
    R = 1 # reflectivity [%] 1 if IAM is IAM_tow(hoy)
    Qin = Ib * R * Trans * IAM * A_helio * Effopt/100 # Eq. 17  in McGovern12
    
    epsilon = 1 # the emissivity of the receiver’s material https://en.wikipedia.org/wiki/Emissivity
    sigma = 5.67 * 1e-8 # Stefan – Boltzman constant [W/m2K4]
    Trec = 565 # the working fluid temperature in the receiver [oC] 565 oC 838K
    Ta = 15 # ambient temperature close to the receiver [oC] 15 oC 288K
    Tin = 290 # working fluid inlet temperature to the receiver [oC] 290 oC 563K
    alpha = 1 # absorptivity of the receiver
    #Ar = 99.3 # receiver area [m2] pp.44 in Pacheco
    
    Qrad = epsilon * sigma * Ar * (CtoK(Trec)**4-CtoK(Ta)**4)
    hconv = CtoK(Trec)/60 + 5/3 # [W/m2K] convection coefficient Eq. 20 in McGovern12
    '''D.L. Siebers, J.S. Kraabel, Estimating convective energy losses from solar central 
    receivers, Sandia National Lab. (SNL-CA), Livermore, CA (United States), 1984. 
    https://doi.org/10.2172/6906848.'''
    Qconv = hconv * Ar * (CtoK(Trec) - CtoK(Ta))
    
    Qnet = alpha * Qin - Qrad - Qconv # Eq. 8 in McGovern12
    
    nR = 0.412 # SAM default

    if Qnet.all() <= 0: #<<<<<<<<<<<<<<<<<<< check with Qin
        P = 0
    else:
        P = Qnet * nR * nG
    return P/1e6 # convert W to MW


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
            Ctow (_type_): the ratio of area of heliostat to solar tower receiver 
        """        
        if not data.shape == (8760,):
            raise ValueError("Wrong time series")
        self.data = data
        self.A_helio = A_helio
        self.Ctow = Ctow

    def hour_power_arr(self):
        return np.vstack((self.hoy, self.data))
    
    def data4surf(self):
        """stacks the data in days and hours 

        Returns:
            _type_: _description_
        """        
        return np.vstack(self.data).reshape((365,24))
    
    @property
    def PowerMax_MW(self)->float:
        """returns the maximum power of the Tower

        Returns:
            _type_: _description_
        """        
        return np.amax(self.data) 

    @property
    def Energy_MWh(self)->float:
        """retunrs the total energy yield of the solar tower. 

        Returns:
           float : the total power per year in [MWh]
        """        
        return integrate.trapz(self.data).round(2)
    @property
    def CF(self)->float: 
        """Capacity factor?

        Returns:
            float: _description_
        """        
        return self.Energy_MWh / (8760 * self.PowerMax_MW)
    
    def as_df(self)->float: 
        return pd.DataFrame({'Power_MW':self.data}, index= self.hoy)

    def summary_tower_data(self):    
        tow_data = np.vstack((self.A_helio,self.Ctow,self.PowerMax_MW,self.Energy_MWh,self.CF)) # vertical stack
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
        self._hourly_results = OutputContainer(data = solarII(Ib=Ib,Trans=transmittance, IAM=IAM_tow(hoy),
                A_helio=self.A_helio_m2,Ar=self.Ar_m2),
            A_helio=self.A_helio_m2, Ctow=self.Ctow)
        return self._hourly_results


#%% Incidence angle methods

def IAM_tow(hoy:np.array=sgh.HOYS_DEFAULT)->float : 
    """Incidence angle modifier of Tower (azimuth)

    for explanation see: http://www.solarpanelsplus.com/solar-tracking/
    
    # polynomial fit, see file IAM.py for data
    
    Args:
        hoy (np.array): hour of year

    Returns:
        float : Incidence angle modifier of Tower in rad
    """    
    return 1.66741484e-1 + 1.41517577e-2 * deg(sgh.z(hoy)) - 9.51787164e-5 * deg(sgh.z((hoy)))**2
    
def IAM_tow2(hoy:np.array=sgh.HOYS_DEFAULT) ->float : # polynomial fit, see file IAM.py for data
    """Incidence angle modifier of Tower - elevation

    Args:
        hoy (np.array): hour of year

    Returns:
        float : Incidence angle modifier of Tower - elevation in rad
    """    
    return 1.66741484e-1 + 1.41517577e-2 * deg(sgh.ele(hoy)) - 9.51787164e-5 * deg(sgh.ele((hoy)))**2



def IAM_tro(hoy:np.array=sgh.HOYS_DEFAULT): 
    """Incidence angle modifier of parabolic trough - equation1
    
    G.A. Salazar, N. Fraidenraich, C.A.A. de Oliveira, O. de Castro Vilela, M. Hongn, J.M. Gordon, 
    Analytic modeling of parabolic trough solar thermal power plants, Energy. 138 (2017) 1148–1156. 
    https://doi.org/10.1016/j.energy.2017.07.110.

    #TODO there are 4 different function for IAM_tro. They need to be consolidated in a single one and selected as an option.

    # thetai in radians

    Args:
        hoy (np.array): hour of year

    Returns:
        _type_: _description_
    """
    #TODO needs rad despite thetai(hoy) already in rad???
    return cos(rad(thetai(hoy))) + 0.02012 * thetai(hoy) - 0.01030 * thetai(hoy)**2 

def IAM_tro2(hoy:np.array=sgh.HOYS_DEFAULT):
    '''N. Fraidenraich, C. Oliveira, A.F. Vieira da Cunha, J.M. Gordon, O.C. Vilela, 
    Analytical modeling of direct steam generation solar power plants, Solar Energy. 98 (2013) 511–522. 
    https://doi.org/10.1016/j.solener.2013.09.037.
    citing M.G. B, E.L. F, R.O. A, A.E. A, W.S. C, A.S. C, E.Z. E, P.N. B, 
    EUROTROUGH- Parabolic Trough Collector Developed for Cost Efficient Solar Power Generation, in: n.d.'''
    #return 1 - 0.00044 * thetai(hoy) / cos(thetai(hoy)) - 0.00003 * (thetai(hoy))**2 / cos(thetai(hoy)) # needs rad despite thetai(hoy) already in rad???
    '''De Luca15
    citing M.J. Montes, A. Abánades, J.M. Martínez-Val, M. Valdés, 
    Solar multiple optimization for a solar-only thermal power plant, using oil as heat transfer fluid in the parabolic trough collectors, 
    Solar Energy. 83 (2009) 2165–2176. https://doi.org/10.1016/j.solener.2009.08.010.
    citing B, M.G.; F, E.L.; A, R.O.; A, A.E.; C, W.S.; C, A.S.; E, E.Z.; B, P.N. 
    EUROTROUGH- Parabolic Trough Collector Developed for Cost Efficient Solar Power Generation.
    '''
    return cos(rad(thetai(hoy))) - 5.25097e-4 * thetai(hoy) - 2.859621e-5 * thetai(hoy)**2

def IAM_tro3(hoy:np.array=sgh.HOYS_DEFAULT):
    '''(Dudley, 1994)'''
    return cos(rad(thetai(hoy))) - 0.0003512 * thetai(hoy) - 0.00003137 * (thetai(hoy))**2# thetai in degrees
    '''pp.26 A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).'''
    #return cos(rad(thetai(hoy))) + 0.000884 * thetai(hoy) - 0.00005369 * (thetai(hoy))**2 # needs rad despite thetai(hoy) already in rad???

def IAM_tro4(hoy,foc_len,area,L):
    '''eq.1 in H.J. Mosleh, R. Ahmadi, Linear parabolic trough solar power plant assisted with latent thermal energy storage system: 
    A dynamic simulation, Applied Thermal Engineering. 161 (2019) 114204. https://doi.org/10.1016/j.applthermaleng.2019.114204.
    J.A. Duffie, W.A. Beckman, Solar Engineering of Thermal Processes, Wiley, 1991.'''
    return 1 - foc_len / L *(1 + area**2 / 48 * foc_len**2) * np.tan(rad(thetai(hoy))) # needs rad despite thetai(hoy) already in rad???


def costhetai(hoy:np.array=sgh.HOYS_DEFAULT): 
    """Parabolic trough cosine function 

    The incidence angle for a plane rotated about a horizontal north-south axis with continuous east
    west tracking to minimize the angle of incidence
    Duffie and Beckman, 1991, pp.24 in  A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).

    # equal to costhetai_NS()
    Args:
        hoy (int):
    
    Returns:
        _type_: _description_
    """    

    return np.sqrt(cos(sgh.z(hoy))**2+np.cos(sgh.d(hoy))**2 * sin(rad(sgh.W(hoy)))**2)

def costhetai_EW(hoy:np.array=sgh.HOYS_DEFAULT):
    """Parabolic trough cosine function in East West orientation

    Gaul, H.; Rabl, A. Incidence-Angle Modifier and Average Optical Efficiency of Parabolic Trough Collectors. 
    Journal of Solar Energy Engineering 1980, 102, 16–21, doi:10.1115/1.3266115.

    Args:
        hoy (int): _description_

    Returns:
        _type_: _description_
    """    
    return np.cos(sgh.d(hoy)) * (cos(rad(W(hoy)))**2 + np.tan(sgh.d(hoy)**2))**0.5

def costhetai_NS(hoy:np.array=sgh.HOYS_DEFAULT):
    """Parabolic trough cosine function in North-South orientation

    Gaul, H.; Rabl, A. Incidence-Angle Modifier and Average Optical Efficiency of Parabolic Trough Collectors. 
    Journal of Solar Energy Engineering 1980, 102, 16–21, doi:10.1115/1.3266115.

    Args:
        hoy (np.array): hour of year

    Returns:
        _type_: _description_
    """    
    return np.cos(sgh.d(hoy)) * (np.sin(rad(sgh.W(hoy)))**2 + (np.cos(rad(sgh.lat)) * np.cos(rad(sgh.W(hoy))) + np.tan(sgh.d(hoy)) * np.sin(rad(sgh.lat)))**2)**0.5


#%% ===========================================================================

def theta_transversal(hoy:np.array=sgh.HOYS_DEFAULT)->float : 
    """Parabolic Trough theta  transversal incidence angle

    Buscemi, A.; Panno, D.; Ciulla, G.; Beccali, M.; Lo Brano, V. 
    Concrete Thermal Energy Storage for Linear Fresnel Collectors: 
    Exploiting the South Mediterranean’s Solar Potential for Agri-Food Processes. 
    Energy Conversion and Management 2018, 166, 719–734, doi:10.1016/j.enconman.2018.04.075.
    
    #TODO  not tested

    Args:
        hoy (np.array): hour of year 

    Returns:
        float: theta  transversal incidence angle
    """    

    return np.arctan(sin(rad(azim(hoy))) * tan(rad(z(hoy))))

def theta_i(hoy:np.array=sgh.HOYS_DEFAULT)->float: 
    """Parabolic Trough longitudinal incidence angle

    Buscemi, A.; Panno, D.; Ciulla, G.; Beccali, M.; Lo Brano, V. 
    Concrete Thermal Energy Storage for Linear Fresnel Collectors: 
    Exploiting the South Mediterranean’s Solar Potential for Agri-Food Processes. 
    Energy Conversion and Management 2018, 166, 719–734, doi:10.1016/j.enconman.2018.04.075.
    
    #TODO  not tested

    Args:
        hoy (np.array): hour of year 

    Returns:
        float: not tested
    """    
    return arctan(cos(rad(azim(hoy))) * tan(rad(z(hoy)))* cos(theta_transversal()))



'''Morin, G.; Dersch, J.; Platzer, W.; Eck, M.; Häberle, A. 
Comparison of Linear Fresnel and Parabolic Trough Collector Power Plants. 
Solar Energy 2012, 86, 1–12, doi:10.1016/j.solener.2011.06.020.'''
# not tested
def thetai_transversal(hoy:np.array=sgh.HOYS_DEFAULT): 
    return np.arctan(abs(sin(azim(hoy)))/tan(ele(hoy)))
def thetai_longtitudinal(hoy:np.array=sgh.HOYS_DEFAULT): 
    return np.arcsin(cos(azim(hoy))*cos(ele(hoy)))

#%% ==================================== solar collector dimenstions
def Ac(Wc:float, L:float, N:int): 
    """returns the collector area. 

    Args:
        Wc (float): width of solar collector in m^2
        L (float): length of solar collector in m^2
        N (int): quantity of solar collectors. 

    Returns:
        _type_: total collector area.
    """    
    return Wc * L * N # update Wc for geometry, see A. Rabl, Comparison of solar concentrators, Solar Energy. 18 (1976) 93–111.

def Ar(Wr, L, N):
    """returns the receiver area. 

    Args:
        Wc (float): width of solar receiver in [m]
        L (float): length of solar receiver in [m]
        N (int): quantity of solar receivers. 

    Returns:
        _type_: total collector area
    """    
    return Wr * L * N

def Cg_tro(Wc, Wr, L, N):
    """Geometrical Concentration 

    Args:
        Wc (float): width of collector in  [m]
        Wr (float): width receiver in  [m]
        L (float): length in  [m]
        N (int): quantity of units []

    Returns:
        _type_: geometrical concentration of parabolic troughs [dimensionless]
    """

    C = Ac(Wc, L, N) / Ar(Wr, L, N)
    return C

def shade_function(Ws,Wc, hoy:np.array=sgh.HOYS_DEFAULT):
    """_summary_

    Args:
        Ws (_type_): #TODO  solar angle
        Wc (_type_): _description_
        hoy (_type_): _description_

    Returns:
        _type_: _description_
    """    
    ''' Stuetzle (2002) pp.29 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).
    N. Fraidenraich, C. Oliveira, A.F. Vieira da Cunha, J.M. Gordon, O.C. Vilela, 
    Analytical modeling of direct steam generation solar power plants, Solar Energy. 98 (2013) 511–522. 
    https://doi.org/10.1016/j.solener.2013.09.037.'''
    return abs(Ws * cos(z(hoy)) / (Wc * cos(thetai(hoy))))

def end_loss(f,L,N, hoy:np.array=sgh.HOYS_DEFAULT):
    '''Lippke, 1995 in pp.31 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).'''
    return (1 - (f * tan(thetai(hoy)) / L)) * N

# def loss_regr(input_dict):
#     '''pp.36-42 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
#     University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).
#     '''
#     equation = a0 + a1 * T + a2 * T**2 + a3 * T**3 + Ib2(alt) * (b0 + b1 * T**2)
#     # Use dict flag to get {variable: value} output, not anonymous [value]
#     #solution = solve(equation.subs(input_dict), dict=True)
#     return equation

#%% total power of parabolic system
def di_sst(Ib,costhetai,IAM,Tr, Wc, Wr, Ws, L, N):
    """Calculates the total power of the parabolic system

    R.K. McGovern, W.J. Smith, Optimal concentration and temperatures of solar thermal power plants,
    Energy Conversion and Management. 60 (2012) 226–232.
    E. Zarza, L. Valenzuela, J. León, K. Hennecke, M. Eck, H.-D. Weyers, M. Eickhoff, 
    Direct steam generation in parabolic troughs: Final results and conclusions of the DISS project, 
    Energy. 29 (2004) 635–644. https://doi.org/10.1016/S0360-5442(03)00172-5.


    Args:
        Ib (_type_): direct irradiance
        costhetai (_type_): cosine function [rad]
        IAM (_type_): incidence angle modifier
        Tr (float): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04 
        Wc (_type_): width of collectors in [m]
        Wr (_type_): width of receiver in [m]
        Ws (_type_): width of spacing between collectors in [m]
        L (_type_): length of units [m]
        N (_type_): number of units

    Returns:
        _type_: power in [MW].
    """    
    Effopt = 75 # [%] Optical efficiency of collector 74% INDITEP in pp.7 Fraidenraich13, pp.4 Zarza06
    
    Qin = Ib * costhetai* IAM * Ac(Wc, L, N) * Effopt/100 # Eq. 4  in McGovern12
    
    epsilon = 1 # the emissivity of the receiver’s material https://en.wikipedia.org/wiki/Emissivity
    sigma = 5.67 * 1e-8 # [W/m2K4] Stefan – Boltzman constant
    #Tr = 350 # [oC] the working fluid temperature in the receiver DISS pp.3,7 in Zarza04
    Ta = 15 # [oC] ambient temperature close to the receiver 15 oC 288K
    Tin = 200 # [oC] working fluid inlet temperature to the receiver DISS pp.3,7 in Zarza04
    alpha = 1 # absorptivity of the receiver

    Qrad = epsilon * sigma * Ar(Wr, L, N) * (CtoK(Tr)**4-CtoK(Ta)**4) # check model from Broesamle
    Qconv = 0
    Qnet = alpha * Qin - Qrad - Qconv # Eq. 8 in McGovern12
    
    # Turbine
    Tcond = 30 # [oC] condenser temperature
    Ts = 565 # [oC] steam temperature
    n_Carnot = 1 - (CtoK(Tcond)/CtoK(Ts)) # Eq. 15  in McGovern12

    nR = 0.375 # isentropic efficiency Salazar17
    if Qnet.all() <= 0:
        P = 0
    else:
        P = Qnet * nR * nG
    return P/1e6 # convert W to MW

def CSCUL(hoys): # CSC heat loss ex.4.2
    """Concentrated Solar Collector heat loss ex.4.2

    Returns:
        _type_: _description_
    """    
    Dr = 0.07 # the receiver’s outer diameter [m]
    Tr = 573 # the working fluid temperature in the receiver [K]
    epsilonr = 0.25 # the emissivity of the receiver’s material
    Dco = 0.1 # the outer diameter of the receiver’s cover [m]
    t = 5*1e-3 # the receiver’s thickness [m]
    u = 4 # wind velocity [m/s]
    Ta = 288 # ambient temperature close to the receiver [K]
    Tsky = 278 # the ambient temperature far from the receiver [K]
    Ucond = 0.022 # the thermal transmittance factor for the heat transfer through conductivity [W/mK]
    epsilonc = 0.90 # the emissivity of the cover’s material
    kc = 1.45 # the receiver’s thermal conductivity factor [W/mK]
    v = 1.456 * 1e-5 # kinematic viscosity of air with ambient temperature 15 οC [m2/s]
    Tco = 306.36 # the temperature of the cover’s outer surface [K]
    sigma = 5.67 * 1e-8 # Stefan – Boltzman constant [W/m2K4]
    keff = 0 # the vacuum’s thermal conductivity factor [W/mK]
    L = 25 # the length of the receiver [m]
    Dci = Dco - (2 * t)
    Re = u * Dco / v
    
    Nu = 0.30 * Re**0.60
    
    hw = Nu * Ucond / Dco
    
    Qloss = pi * Dco * hw * (Tco - Ta) 
    + epsilonc * pi * Dco * sigma * (Tco**4 - Tsky**4) * L
    Tci = Qloss / (2 * pi * kc * L) * ln(Dco/Dci) + Tco
    Qloss2 = (2 * pi * keff * (CtoK(Tr) - Tci) / ln(Dci/Dr) + 
    (pi * Dr * sigma * (CtoK(Tr)**4 - Tci**4)) / (1/epsilonr+((1-epsilonc)*Dr/Dci)/epsilonc)) * L
    for x in hoys:
        UL = Qloss / (pi * Dr * L * (CtoK(Tr) - Ta))
    return UL # [W/m2K]

def CSCP(Tfi, hoy:np.array= sgh.HOYS_DEFAULT, fname:str="example_data/tmy_35.015_25.755_2005_2020.csv"): 
    # CSC thermal power ex.4.3
    """Concentrated Solar collect

    Args:
        Tfi (_type_): _description_
        hoy (np.array, optional): _description_. Defaults to sgh.HOYS_DEFAULT.
        fname (str, optional): _description_. Defaults to "example_data/tmy_35.015_25.755_2005_2020.csv".

    Returns:
        _type_: _description_
    """    
    '''see also Dikmen, E., Ayaz, M., Ezen, H.H., Küçüksille, E.U., Şahin, A.Ş., 2014. 
    Estimation and optimization of thermal performance of evacuated tube solar collector system. 
    Heat Mass Transfer 50, 711–719. https://doi.org/10.1007/s00231-013-1282-0
    '''
    
    try:
        pvgis = pd.read_csv(fname, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python')
    except:
        #TODO this is an exception until T = CSCP(Tr) is removed from this file (use tests)
        pvgis = pd.read_csv("examples/"+fname, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python')

    Ib = pvgis.loc[:,'Gb(n)']
    #S = 550 # incident solar radiation [W/m2]
    #UL = 4.50 # heat losses total thermal transmittance factor [W/m2Κ]
    Wc = 5.76 # width of the concentrating collector [m]
    L = 12 # length of the concentrating collector [m]
    cp = 2500 # fluid's specific heat capacity [J/kgK]
    m = 8 # flow rate [kg/s] Mosleh19, pp.9
    #Tfi = Tr # fluid's inlet temperature [C]
    hfi = 400 # thermal convection factor for the heat transfer from the receiver to the working fluid [W/m2K]
    k = 16 # stainless steel thermal conductivity factor [W/mK]
    t = 5 * 1e-3 # receiver's tube thickness [m]
    Ta = 15 # ambient temperature close to the receiver [C]
    Dro = 0.07 # receiver’s outer diameter [m]
    Dco = 0.1 # outer diameter of the receiver’s cover [m]
    Ti = 270 #  average temperature of the receiver [C]
    
    Ar = pi * Dro * L # receiver’s inner area [m2]
    Aa = (Wc - Dco) * L # collector’s effective area [m2]
    Dri = Dro - (2 * t) #  receiver’s inner diameter [m]
    F = 1/CSCUL(hoy) / (1/CSCUL(hoy) + Dro/(hfi*Dri) + Dro * ln(Dro/Dri)/(2*k))
    FR = m * cp * (1 - np.exp(-Ar * CSCUL(hoy) * F / (m * cp))) / (Ar * CSCUL(hoy)) # heat removal factor - review precision <<<
    Qu = Aa * FR * (Ib - Ar * CSCUL(hoy) * (Ti - Ta)/Aa) # final thermal power production [W]
    Tfo = Tfi + (Qu / (m * cp)) # fluid’s outlet temperature [C]
    Tro =  Tfi + (Qu * ((1/(pi * Dri * L * hfi)) + ln(Dro/Dri)/(2*pi*k*L))) # receiver’s outer surface temperature [C]
    DT = Tro - Tfi
    return Tfo#Qu/1000 # convert W to kW




# %%
