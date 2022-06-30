#%%
'''Concentrating Solar Power plants                 ALT + SHIFT +0 to unfold levels
                                                    ALT + 0 to fold levels'''
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
# from pylab import *

# from SolarGeometry_hoy import *
import SolarGeometry_hoy as sgh
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
    
def solarII(Ib,Trans,IAM,A_helio:float,Ar:float)->float:
    """Calculates the power of the solar tower
 
    R.K. McGovern, W.J. Smith, Optimal concentration and temperatures of solar thermal power plants,
    Energy Conversion and Management. 60 (2012) 226–232.
    J.E. Pacheco, R.W. Bradshaw, D.B. Dawson, W.D. la Rosa, R. Gilbert, S.H. Goods, P. Jacobs, M.J. Hale, S.A. Jones, G.J. Kolb, M.R. Prairie, H.E. Reilly, S.K. Showalter, L.L. VANT-HULL, 
    Final Test and Evaluation Results from the Solar Two Project, n.d. https://core.ac.uk/reader/193342950 (accessed September 8, 2020).
 


    Args:
        Ib (_type_): direct irradiance
        Trans (_type_): transmissivity 
        IAM (_type_): incidence angle modifier
        A_helio (float): heliostat area in m^2
        Ar (float): receiver area in m^2

    Returns:
        float: power in MW
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

def IAM_tow(hoy:np.array=sgh.HOYS_DEFAULT)->float : 
    """Incidence angle modifier of Tower (azimuth)
    
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
        Tr (_type_): Tranmissivity 
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

def CSCP(Tfi, hoy:np.array= sgh.HOYS_DEFAULT, fname:str="example_data/tmy_35.015_25.755_2005_2020.csv"): # CSC thermal power ex.4.3
    '''see also Dikmen, E., Ayaz, M., Ezen, H.H., Küçüksille, E.U., Şahin, A.Ş., 2014. 
    Estimation and optimization of thermal performance of evacuated tube solar collector system. 
    Heat Mass Transfer 50, 711–719. https://doi.org/10.1007/s00231-013-1282-0
    '''
    
    pvgis = pd.read_csv(fname, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python')
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

def heatloss_Mertins(d,epsilon,dt):
    '''eq.20 in [1]M.J. Montes, R. Barbero, R. Abbas, A. Rovira, 
    Performance model and thermal comparison of different alternatives for the Fresnel single-tube receiver, 
    Applied Thermal Engineering. 104 (2016) 162–175. https://doi.org/10.1016/j.applthermaleng.2016.05.015.
    '''
    return (d / 0.219) * ((1.945 - 0.2428 * epsilon/0.08) * dt + (0.001226 + 0.004568 * epsilon/0.08) * dt**2)

def pipe_loss():
    '''Price, 2005 in pp.42 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).
    '''
    DT = (sgh.Tout + sgh.Tin)/2 - sgh.Ta # [oC]
    # 0.017 * DT - 1.683 * 1e-4 * DT**2 + 6.78 * 1e-7 * DT**3 ref?
    return 0.01693 * DT - 0.0001683 * DT**2 + 6.78 * 1e-7 * DT**3

#PCM data
pcm_data = [[2200,2257,2110,2044,2380],
            [212,174,226,149.7,280],
            [282,308,333,380,250,802],
            [NaN,0.5,0.5,0.5,0.52,5.0],
            [1.733,1.588,1.240,NaN,NaN,NaN],
            [2.553,1.650,1.341,NaN,NaN,NaN],
            [0.2,0.2,0.3,1.0,NaN,0.15]] #NaNO2 assumption
pcm = pd.DataFrame(pcm_data,
    columns = ['NaNO2','NaNO3','KNO3','KOH','H250','NaCl'])
pcm_rho = pcm['NaNO3'].loc[0]
pcm_latent_heat = pcm['NaNO3'].loc[1]#*0.27778 # [kJ/kgK to Wh/kgK] https://www.cactus2000.de/uk/unit/masscp1.php
pcm_melting_point = pcm['NaNO3'].loc[2] # [oC]
pcm_thermal_conductivity_coeff = pcm['NaNO3'].loc[3]
pcm_solid_cp = pcm['NaNO3'].loc[4] # [kJ/kgK]
pcm_liquid_cp = pcm['NaNO3'].loc[5] # [kJ/kgK]
pcm_cost = pcm['NaNO3'].loc[6]

def nano3_kno3_cp():
    '''G.J. Janz, C.B. Allen, N.P. Bansal, R.M. Murphy, R.P.T. Tomkins, 
    Physical Properties Data Compilations Relevant to Energy Storage. II. 
    Molten Salts: Data on Single and Multi-Component Salt Systems, 
    NATIONAL STANDARD REFERENCE DATA SYSTEM, 1979. 
    https://apps.dtic.mil/sti/citations/ADD095217 (accessed April 21, 2021).
    '''
    # pp.396-406(405-415) in Janz79 Tmin=222C
    nano3_kno3 = [[0,25,40,50,60,75,87,100],
                  [3.600,3.382,3.290,3.195,3.100,2.802,2.685,2.300],
                  [310,275,243,227,230,260,300,337]] # [mol%kno3 vs DH kcal/mol vs Tm oC]
    mol = [0,25,40,50,60,75,87,100]
    dh = [3.600,3.382,3.290,3.195,3.100,2.802,2.685,2.300]
    T = [310,275,243,227,230,260,300,337]
    plt.plot(mol,dh)
    plt.plot(mol,T)
    for T in np.arange(510,770,10):
        nano3_kno3_cp = 53.44 - 2.638 * 10**-2 * T # [T in K, cp in cal/K mol]
        plt.scatter(T, nano3_kno3_cp)

def delta_h(a,b,c,x_a,x_b):
    '''Coscia, K.; Elliott, T.; Mohapatra, S.; Oztekin, A.; Neti, S. 
    Binary and Ternary Nitrate Solar Heat Transfer Fluids. 
    Journal of Solar Energy Engineering 2013, 135, 
    doi:10.1115/1.4023026.
    '''
    # [cal/mole]
    nano3_kno3_l = (-1707, -284, 0)
    nano3_kno3_s = (6276, 0, 0)
    lino3_nano3_l = (-1941, -2928, 0)
    lino3_nano3_s = (9204, 3347, 0)
    lino3_kno3_l = (-9183, -364, -1937)
    lino3_kno3_s = (10460, 4184, 0)
    return x_a * x_b * (a + b * x_a + c * x_a * x_b)

def delta_s(x_a):
    '''Coscia, K.; Elliott, T.; Mohapatra, S.; Oztekin, A.; Neti, S. 
    Binary and Ternary Nitrate Solar Heat Transfer Fluids. 
    Journal of Solar Energy Engineering 2013, 135, 
    doi:10.1115/1.4023026.
    '''
    R = 8.31446261815324 # [J⋅K−1⋅mol−1]
    return R * ln(x_a)

def mass_mix(x_a,x_b,m_a,m_b):
    m_nano3 = 84.9947 #g/mol
    m_kno3 = 101.1032 #g/mol
    m_lino3 = 68.946 #g/mol
    return x_a * m_a + x_b * m_b

def cp_mix():
    '''
    B. D’Aguanno, M. Karthik, A.N. Grace, A. Floris, 
    Thermostatic properties of nitrate molten salts and their solar and eutectic mixtures, 
    Sci Rep. 8 (2018) 10485. https://doi.org/10.1038/s41598-018-28641-1.
    '''
    c1 = 1.81 #NaNO3
    c2 = 1.52 # KNO3
    for x1 in np.arange(0,1,0.1):
        x2 = 1 - x1
        rmix = x1 * c1 + x2 * c2
        plt.scatter(x1, rmix)
    '''
    M. Liu, W. Saman, F. Bruno, 
    Review on storage materials and thermal performance enhancement techniques for high temperature 
    phase change thermal storage systems, Renewable and Sustainable Energy Reviews. 16 (2012) 2118–2132. 
    https://doi.org/10.1016/j.rser.2012.01.020.
    '''
    pcm_solid_cp1 = pcm['NaNO2'].loc[4] # [kJ/kgK]
    pcm_liquid_cp1 = pcm['NaNO2'].loc[5] # [kJ/kgK]
    pcm_solid_cp2 = pcm['NaNO3'].loc[4] # [kJ/kgK]
    pcm_liquid_cp2 = pcm['NaNO3'].loc[5] # [kJ/kgK]
    pcm_solid_cp3 = pcm['KNO3'].loc[4] # [kJ/kgK]
    pcm_liquid_cp3 = pcm['KNO3'].loc[5] # [kJ/kgK]

    for x1 in np.arange(0,1,0.1):
        x2 = 1 - x1
        #1.15 scale to 
        pcm_solid_cpmix1 = x1 * pcm_solid_cp1+ x2 * pcm_solid_cp2
        pcm_liquid_cpmix1 = x1 * pcm_liquid_cp1 + x2 * pcm_liquid_cp2
        plt.scatter(x1, pcm_solid_cpmix1, c='b')
        plt.scatter(x1, pcm_liquid_cpmix1, c='g')
    
    for x1 in np.arange(0,1,0.1):
        x2 = 1 - x1
        #1.15 scale to [D’Aguanno18] experiments
        pcm_solid_cpmix2 = x1 * pcm_solid_cp2*1.15 + x2 * pcm_solid_cp3*1.15 
        pcm_liquid_cpmix2 = x1 * pcm_liquid_cp2*1.15 + x2 * pcm_liquid_cp3*1.15
        plt.scatter(x1, pcm_solid_cpmix2, c='c')
        plt.scatter(x1, pcm_liquid_cpmix2, c='r')
    return

# define materials
Vpcm = 200 # [m3]
#mpcm = pcm_rho * Vpcm # kg/m3 * m3
mpcm = 13.7 # [kg/s] Mahfuz14 pp.5
hp = pcm_latent_heat
Ap = (pow(Vpcm,1/3))**2 # [m2] area of storage

cs = pcm_solid_cp #* 2.7778e-4 # [J/kgK to Wh/kgK] https://www.cactus2000.de/uk/unit/masscp1.php
ce = pcm_latent_heat #* 2.7778e-4 # [J/kgK to Wh/kgK]
cl = pcm_liquid_cp #* 2.7778e-4 # [J/kgK to Wh/kgK]

Tr = pcm_melting_point + 10 # [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04
Tmelt = pcm_melting_point
T2 = pcm_melting_point + 5 #Mahfuz14 pp.5
#T2 = pcm_temp()
T = CSCP(Tr)

def pcm_temp():
    '''https://stackabuse.com/solving-systems-of-linear-equations-with-pythons-numpy/
    Y.-Q. Li, Y.-L. He, Z.-F. Wang, C. Xu, W. Wang, Exergy analysis of two phase change materials storage system 
    for solar thermal power with finite-time thermodynamics, Renewable Energy. 39 (2012) 447–454. 
    https://doi.org/10.1016/j.renene.2011.08.026.
    H.J. Mosleh, R. Ahmadi, Linear parabolic trough solar power plant assisted with latent thermal energy storage system: 
    A dynamic simulation, Applied Thermal Engineering. 161 (2019) 114204. https://doi.org/10.1016/j.applthermaleng.2019.114204.
    '''

    Tin = Tr
    Nc = hp * Ap / mpcm * cpcm
    T2 = Tmelt + (T - Tmelt) * np.exp(-Nc)
    return T2 #Estore

def phase_change(T):
    '''L. Solomon, A. Oztekin, Exergy analysis of cascaded encapsulated phase change material—
    High-temperature thermal energy storage systems, Journal of Energy Storage. 8 (2016) 12–26. 
    https://doi.org/10.1016/j.est.2016.09.004.
    M. Thonon, G. Fraisse, L. Zalewski, M. Pailha, 
    Towards a better analytical modelling of the thermodynamic behaviour of phase change materials, 
    Journal of Energy Storage. 32 (2020) 101826. https://doi.org/10.1016/j.est.2020.101826.
    '''
    #exp_rat = exp() / Tl-Ts/2 * sqrt(pi)
    #ceff_s = cs + pcm_latent_heat * exp_rat
    #ceff_s = cl + pcm_latent_heat * exp_rat
    if T[0] < Tmelt:
        gamma_fraction = 0
        #Qstor = mpcm*cs*(T - Tmelt)
        Qstor = mpcm*cs*(T - Tmelt)
    elif T[0] >= Tmelt and T[0] <= T2:
        gamma_fraction = (T - Tmelt) / (T2 - Tmelt)
        #Qstor = mpcm*cs*(Tmelt - T) + mpcm*ce*(T - Tmelt)
        Qstor = mpcm*cs*(T - Tmelt) + mpcm*gamma_fraction*pcm_latent_heat
    elif T[0] > T2:
        gamma_fraction = 1
        #Qstor = mpcm*cs*(Tmelt - T) + mpcm*ce*(T - Tmelt) + mpcm*cl*(T - T2)
        Qstor = mpcm*cs*(T - Tmelt) + mpcm*gamma_fraction*pcm_latent_heat + mpcm*cl*(T - Tmelt)
    #plt.plot(T, gamma_fraction)
    return Qstor/1e3 # [kW to MW] gamma_fraction,{'gamma_fraction': gamma_fraction, 'Qstor': Qstor}

# %%
