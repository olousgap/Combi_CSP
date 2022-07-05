# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis, G. Arnaoutakis
    @Date: 2022/07/02
    @Credit: original functions from G. Arnaoutakis
"""
#%%
from multiprocessing.sharedctypes import Value
import numpy as np
import pandas as pd

from CombiCSP import OutputContainer, CtoK, HOYS_DEFAULT
import CombiCSP.SolarGeometry as sgh
import CombiCSP.CSP as cspC

from CombiCSP.SolarGeometry import thetai, W
from CombiCSP.solar_system_location import SolarSystemLocation, d


class SolarTroughCalcs():
    def __init__(self,
        foc_len = 0.88 
        ,N = 1800 
        ,L = 25 
        ,Ws = 18 
        ,Wr = 0.07 
        ,Wc = 5.76
        , slobj:SolarSystemLocation =  None
        ):
        """_summary_

        Args:
            foc_len (float): [m] focal length CSPP T.1 in Mosleh19. Defaults to 0.88.
            N (int):  # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine. Defaults to 1800.
            L (int): [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine  . Defaults to 25.
            Ws (int): [m] width between rows 18 INDITEP in pp.6 Fraidenraich13,  pp.5 Zarza06. Defaults to 18.
            Wr (float): tube outer diameter [m]. Defaults to 0.07.
            Wc (float): collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM. Defaults to 5.76.
        """        
        if isinstance(slobj, SolarSystemLocation):
            self._sl = slobj
        else:
            raise ValueError('Site Location has not been provided')

        self.foc_len = foc_len
        self.N = N
        self.L = L
        self.Ws = Ws
        self.Wr = Wr
        self.Wc = Wc 

    def print_system_summary(self):
        print(f"""
--------  System Parameters
-> Focal Lentth     (foc_len) = {self.foc_len}  
-> Number of units  (   N   ) = {self.N}      
-> Unit length      (   L   ) =  {self.L} 
-> Unit spacing?    (   Ws  ) = {self.Ws} 
-> Receiver width   (   Wr  ) = {self.Wr}
-> Collector width  (   Wc  ) = {self.Wc}
--------  Derived Quantities
-> Collector area   (   Ac  ) = {self.Ac()}
-> Receiver area    (   Ar  ) = {self.Ar()}
-> Concetration f   (   Cg  ) = {self.Cg}
        """)        

    @property
    def area(self):
        """returns the collector area. 

        update Wc for geometry, see A. Rabl, Comparison of solar concentrators, Solar Energy. 18 (1976) 93–111.
        
        Args:
            Wc (float): width of solar collector in m^2
            L (float): length of solar collector in m^2
            N (int): quantity of solar collectors. 

        Returns:
            _type_: total collector area.


        Returns:
            float: The collector area in [m^2]
        """        
        return self.Wc * self.L * self.N
    
    def Ac(self)->float:
        """collector area in m^2

        (Assumption): it is the width*length times the number of units

        Returns:
            float: collector area in m^2
        """        
        return self.Wc * self.L * self.N

    def Ar(self)->float:
        """returns the receiver area. 

        Args:
            Wr (float): width of solar receiver in [m]
            L (float): length of solar receiver in [m]
            N (int): quantity of solar receivers. 

        Returns:
            float: total receiver area in m^2
        """    
        return self.Wr * self.L * self.N
    
    @property
    def Cg(self)->float:
        """Geometrical Concentration 

        Affected by the following properties:
            Wc (float): width of collector in  [m]
            Wr (float): width receiver in  [m]
            L (float): length in  [m]
            N (int): quantity of units []

        Returns:
            float: geometrical concentration of parabolic troughs [dimensionless]
        """

        return self.Ac() / self.Ar()

    def perform_calcs_EW(self, Ib, Tr=318, hoy=HOYS_DEFAULT):
        """Calculation for a solar trough oriented EW for a year per hour 

        Args:
            Ib (pd.Series): beam irradiance
            Tr (float, optional): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04. Defaults to 318.
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.

        Returns:
            OutputContainer: Object that contains the power [MW] for each hour for the trough.
        """ 
        IAM = self.IAM_tro(hoy)
        
        #Parabolic trough cosine function in East West orientation
        #    Gaul, H.; Rabl, A. Incidence-Angle Modifier and Average Optical Efficiency of Parabolic Trough Collectors. 
        #   Journal of Solar Energy Engineering 1980, 102, 16–21, doi:10.1115/1.3266115.
        costhetai_EW_arr =  np.cos( d(hoy)) * (np.cos(np.radians(self._sl.W(hoy)))**2 + np.tan(d(hoy)**2))**0.5
        
        data = self.di_sst(hoy = hoy, Ib=Ib,costhetai= costhetai_EW_arr, Tr=Tr)
        return OutputContainer(data = data, A_helio=self.area, Ctow=self.Cg)

    def costhetai_EW(self, hoy):
        return  np.cos( d(hoy)) * (np.cos(np.radians(self._sl.W(hoy)))**2 + np.tan(d(hoy)**2))**0.5

    def perform_calcs_NS(self, Ib, Tr=318., hoy=HOYS_DEFAULT):
        """Calculation for a solar trough oriented NS for a year per hour 

        Args:
            Ib (pd.Series): beam irradiance
            Tr (float, optional): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04. Defaults to 318.
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.

        Returns:
            OutputContainer: Object that contains the power [MW] for each hour for the trough.
        """        

        lat_rad = self._sl.lat_rad
        #Parabolic trough cosine function in North-South orientation
        #   Gaul, H.; Rabl, A. Incidence-Angle Modifier and Average Optical Efficiency of Parabolic Trough Collectors. 
        #   Journal of Solar Energy Engineering 1980, 102, 16–21, doi:10.1115/1.3266115.
        costhetai_NS_arr = np.cos(d(hoy)) * (np.sin(np.radians(self._sl.W(hoy)))**2 + 
                (np.cos(lat_rad) *  np.cos(np.radians(self._sl.W(hoy))) + np.tan(d(hoy)) * np.sin(lat_rad))**2)**0.5
 
        data = self.di_sst(hoy=hoy, Ib=Ib,costhetai=costhetai_NS_arr,
                      Tr=Tr)
        return OutputContainer(data = data, A_helio=self.area, Ctow=self.Cg)
    
    def costhetai_NS(self, hoy):
        lat_rad = self._sl.lat_rad
        return np.cos(d(hoy)) * (np.sin(np.radians(self._sl.W(hoy)))**2 + 
                (np.cos(lat_rad) *  np.cos(np.radians(self._sl.W(hoy))) + np.tan(d(hoy)) * np.sin(lat_rad))**2)**0.5
 

    def thetai(self, hoy:np.array=HOYS_DEFAULT, inclination=90., azimuths=0.)->np.array: #
        """ Calculates the incidence angle [in radians]

        #TODO check whether there is a dependence between azimuth and NS and EW type of CSP

        Args:
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
            inclination (float, optional): _description_. Defaults to 90.
            azimuths (float, optional): _description_. Defaults to 0.

        Returns:
            np.array: _description_
        """        
        g = np.degrees(self._sl.azim(hoy)) - azimuths # if surface looks due S then azimuths=0
        elev = self._sl.ele(hoy)
        return np.arccos(np.cos(elev) * np.sin(np.radians(inclination)) * np.cos(np.radians(g)) 
            + np.sin(elev) * np.cos(np.radians(inclination)))


    def IAM_tro(self, hoy:np.array=HOYS_DEFAULT): 
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
        return np.cos(np.radians(self.thetai(hoy))) + 0.02012 * self.thetai(hoy) - 0.01030 * self.thetai(hoy)**2 

    def di_sst(self, hoy, Ib, costhetai, Tr, nG:float = 0.97)->pd.Series:
        """Calculates the total power of the parabolic system

        R.K. McGovern, W.J. Smith, Optimal concentration and temperatures of solar thermal power plants,
        Energy Conversion and Management. 60 (2012) 226–232.
        E. Zarza, L. Valenzuela, J. León, K. Hennecke, M. Eck, H.-D. Weyers, M. Eickhoff, 
        Direct steam generation in parabolic troughs: Final results and conclusions of the DISS project, 
        Energy. 29 (2004) 635–644. https://doi.org/10.1016/S0360-5442(03)00172-5.


        Args:
            hoy (np.array): hour of year
            Ib (np.array): hour of year
            costhetai (_type_): cosine function [rad]
            IAM (_type_): incidence angle modifier
            Tr (float): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04 
            nG (float): Generator efficiency [dimensionless]

        Returns:
            _type_: power in [MW].
        """    
        Wc = self.Wc # width of collectors in [m]
        Wr = self.Wr # width of receiver in [m]
        Ws = self.Ws # width of spacing between collectors in [m]
        L = self.L   # length of units [m]
        N = self.N   # number of units
        
        IAM = self.IAM_tro(hoy) # incidence angle modifier 
        
        
        Effopt = 75 # [%] Optical efficiency of collector 74% INDITEP in pp.7 Fraidenraich13, pp.4 Zarza06
        
        Qin = Ib * costhetai* IAM * self.Ac() * Effopt/100 # Eq. 4  in McGovern12
        
        epsilon = 1 # the emissivity of the receiver’s material https://en.wikipedia.org/wiki/Emissivity
        sigma = 5.67 * 1e-8 # [W/m2K4] Stefan – Boltzman constant
        #Tr = 350 # [oC] the working fluid temperature in the receiver DISS pp.3,7 in Zarza04
        Ta = 15 # [oC] ambient temperature close to the receiver 15 oC 288K
        Tin = 200 # [oC] working fluid inlet temperature to the receiver DISS pp.3,7 in Zarza04
        alpha = 1 # absorptivity of the receiver

        Qrad = epsilon * sigma * self.Ar() * (CtoK(Tr)**4-CtoK(Ta)**4) # check model from Broesamle
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

    
    def mutate(self,foc_len = None 
        ,N = None 
        ,L = None 
        ,Ws = None
        ,Wr = None 
        ,Wc = None
        , slobj:SolarSystemLocation =  None):
        foc_len = self.foc_len if foc_len is None else foc_len
        N = self.N if N is None else N
        L = self.L if L is None else L
        Ws = self.Ws if Ws is None else Ws
        Wr = self.Wr if Wr is None else Wr
        Wc = self.Wc if Wc is None else Wc
        slobj = self._sl if slobj is None else slobj
        return SolarTroughCalcs(
            foc_len = foc_len 
            ,N = N 
            ,L = L 
            ,Ws = Ws 
            ,Wr = Wr 
            ,Wc =Wc 
            , slobj=  slobj
                                )
        

#%% Incidence angle methods for troughs

def IAM_tro(hoy:np.array=HOYS_DEFAULT): 
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
    return np.cos(np.radians(thetai(hoy))) + 0.02012 * thetai(hoy) - 0.01030 * thetai(hoy)**2 

def IAM_tro2(hoy:np.array=HOYS_DEFAULT):
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
    return np.cos(np.radians(thetai(hoy))) - 5.25097e-4 * thetai(hoy) - 2.859621e-5 * thetai(hoy)**2

def IAM_tro3(hoy:np.array=HOYS_DEFAULT):
    '''(Dudley, 1994)'''
    return np.cos(np.radians(thetai(hoy))) - 0.0003512 * thetai(hoy) - 0.00003137 * (thetai(hoy))**2# thetai in degrees
    '''pp.26 A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).'''
    #return cos(rad(thetai(hoy))) + 0.000884 * thetai(hoy) - 0.00005369 * (thetai(hoy))**2 # needs rad despite thetai(hoy) already in rad???

def IAM_tro4(hoy,foc_len,area,L):
    '''eq.1 in H.J. Mosleh, R. Ahmadi, Linear parabolic trough solar power plant assisted with latent thermal energy storage system: 
    A dynamic simulation, Applied Thermal Engineering. 161 (2019) 114204. https://doi.org/10.1016/j.applthermaleng.2019.114204.
    J.A. Duffie, W.A. Beckman, Solar Engineering of Thermal Processes, Wiley, 1991.'''
    return 1 - foc_len / L *(1 + area**2 / 48 * foc_len**2) * np.tan(np.rad(thetai(hoy))) # needs rad despite thetai(hoy) already in rad???


def costhetai(hoy:np.array=HOYS_DEFAULT): 
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

    return np.sqrt(np.cos(sgh.z(hoy))**2+np.cos(sgh.d(hoy))**2 * np.sin(np.radians(sgh.W(hoy)))**2)

def costhetai_EW(hoy:np.array=HOYS_DEFAULT):
    """Parabolic trough cosine function in East West orientation

    Gaul, H.; Rabl, A. Incidence-Angle Modifier and Average Optical Efficiency of Parabolic Trough Collectors. 
    Journal of Solar Energy Engineering 1980, 102, 16–21, doi:10.1115/1.3266115.

    Args:
        hoy (int): _description_

    Returns:
        _type_: _description_
    """    
    return np.cos(sgh.d(hoy)) * (np.cos(np.radians(W(hoy)))**2 + np.tan(sgh.d(hoy)**2))**0.5

def costhetai_NS(hoy:np.array=HOYS_DEFAULT):
    """Parabolic trough cosine function in North-South orientation

    Gaul, H.; Rabl, A. Incidence-Angle Modifier and Average Optical Efficiency of Parabolic Trough Collectors. 
    Journal of Solar Energy Engineering 1980, 102, 16–21, doi:10.1115/1.3266115.

    Args:
        hoy (np.array): hour of year

    Returns:
        _type_: _description_
    """    
    return np.cos(sgh.d(hoy)) * (np.sin(np.radians(sgh.W(hoy)))**2 + (np.cos(np.radians(sgh.lat)) * np.cos(np.radians(sgh.W(hoy))) + np.tan(sgh.d(hoy)) * np.sin(np.radians(sgh.lat)))**2)**0.5

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


#%% total power of parabolic system
def di_sst(Ib,costhetai,IAM,Tr, Wc, Wr, Ws, L, N, nG:float = 0.97):
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
    
    Qloss = np.pi * Dco * hw * (Tco - Ta) 
    + epsilonc * np.pi * Dco * sigma * (Tco**4 - Tsky**4) * L
    Tci = Qloss / (2 * np.pi * kc * L) * np.log(Dco/Dci) + Tco
    Qloss2 = (2 * np.pi * keff * (CtoK(Tr) - Tci) / np.log(Dci/Dr) + 
    (np.pi * Dr * sigma * (CtoK(Tr)**4 - Tci**4)) / (1/epsilonr+((1-epsilonc)*Dr/Dci)/epsilonc)) * L
    for x in hoys:
        UL = Qloss / (np.pi * Dr * L * (CtoK(Tr) - Ta))
    return UL # [W/m2K]

def CSCP(Tfi, hoy:np.array= HOYS_DEFAULT, fname:str="example_data/tmy_35.015_25.755_2005_2020.csv"): 
    # CSC thermal power ex.4.3
    """Concentrated Solar collect

    Args:
        Tfi (_type_): _description_
        hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
        fname (str, optional): _description_. Defaults to "example_data/tmy_35.015_25.755_2005_2020.csv".

    Returns:
        _type_: _description_
    """    
    '''see also Dikmen, E., Ayaz, M., Ezen, H.H., Küçüksille, E.U., Şahin, A.Ş., 2014. 
    Estimation and optimization of thermal performance of evacuated tube solar collector system. 
    Heat Mass Transfer 50, 711–719. https://doi.org/10.1007/s00231-013-1282-0
    '''
    
    try:
        pvgis_data = pd.read_csv(fname, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python')
    except:
        #TODO this is an exception until T = CSCP(Tr) is removed from this file (use tests)
        pvgis_data = pd.read_csv("examples/"+fname, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python')

    Ib = pvgis_data.loc[:,'Gb(n)']
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
    
    Ar = np.pi * Dro * L # receiver’s inner area [m2]
    Aa = (Wc - Dco) * L # collector’s effective area [m2]
    Dri = Dro - (2 * t) #  receiver’s inner diameter [m]
    F = 1/CSCUL(hoy) / (1/CSCUL(hoy) + Dro/(hfi*Dri) + Dro * np.log(Dro/Dri)/(2*k))
    FR = m * cp * (1 - np.exp(-Ar * CSCUL(hoy) * F / (m * cp))) / (Ar * CSCUL(hoy)) # heat removal factor - review precision <<<
    Qu = Aa * FR * (Ib - Ar * CSCUL(hoy) * (Ti - Ta)/Aa) # final thermal power production [W]
    Tfo = Tfi + (Qu / (m * cp)) # fluid’s outlet temperature [C]
    Tro =  Tfi + (Qu * ((1/(np.pi * Dri * L * hfi)) + np.log(Dro/Dri)/(2*np.pi*k*L))) # receiver’s outer surface temperature [C]
    DT = Tro - Tfi
    return Tfo#Qu/1000 # convert W to kW

