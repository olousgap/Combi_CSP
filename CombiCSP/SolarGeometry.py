# -*- coding: utf-8 -*- 
"""
    @Author: G.Arnaoutakis
    @Date: 2022/mm/yy
    @Credit: heavily modified by N. Papadakis on 2022/07/22
"""
#%%
'''Geometry calculations as a function of hourly local time '''
import numpy as np
import pandas as pd
# import pvlib
# from pvlib import clearsky, atmosphere, solarposition

from CombiCSP import HOYS_DEFAULT
from CombiCSP.solar_system_location import EoT, d 




#TODO Move this parameters to another location because now the system can only calculate for these long and lat.
dt_gmt = +2 # time difference between Greenwich Mean Time
lat = 35 # Crete
mer = -25 # for Greece check to replace with 15 * dt_gmt
lon = 24 # Crete 35.2401° N, 24.8093° E [east negative, west positive]



#%% 
def tsol(hoy:np.array=HOYS_DEFAULT): # solar time [in decimal hours] introduce if function for east/west<<<<<<<<<
    """returns solar time 

    CALCULATION NEGLECTS DAYLIGHT SAVING ON SUMMER
    https://www.pveducation.org/pvcdrom/properties-of-sunlight/the-suns-position
    hoy + ((lat-mer)/15 + EoT(hoy))/60 # [60 min/h]
    pp.5, Appendix C in D.A. Katsaprakakis, Power Plant Synthesis, CRC Press, 2020.
    https://doi.org/10.1201/b22190

    Args:
        hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.

    Returns:
        _type_: _description_
    """    
    return hoy + (4*(lon-15*dt_gmt) + EoT(hoy))/60 # [60 min/h]


def W(hoy:np.array=HOYS_DEFAULT): # solar hour angle [in degrees]
    '''given than for tsol = 12h it should be ω = 0ο and for the solar time range
    tsol = 0 – 24h the solar hour angle ranges from 0o to ±180'''
    return 15 * (tsol(hoy) - 12) # 360deg/24h = 15deg/h

def ele(hoy:np.array=HOYS_DEFAULT): # solar elevation angle or solar height [in radians]
    return np.arcsin(np.cos(np.radians(lat)) * np.cos(d(hoy)) * np.cos(np.radians(W(hoy)))
    + np.sin(np.radians(lat)) * np.sin(d(hoy)))

def z(hoy:np.array=HOYS_DEFAULT)->np.array:
    """Returns the solar zenith angle in radians

    Args:
        hoy (np.array, optional): hour of year. Defaults to HOYS_DEFAULT.

    Returns:
        np.array: solar zenith angle in radians
    """    
    # solar zenith angle [in radians] https://en.wikipedia.org/wiki/Solar_zenith_angle
    return np.arccos(np.cos(np.radians(lat)) * np.cos(d(hoy)) * np.cos(np.radians(W(hoy))) \
        + np.sin(np.radians(lat)) * np.sin(d(hoy)))

def azim(hoy:np.array=HOYS_DEFAULT)->np.array: 
    """Returns the solar azimuth angle in radians

    Args:
        hoy (np.array, optional): hour of year. Defaults to HOYS_DEFAULT.

    Returns:
        np.array: solar azimuth angle in radians
    """   
    return np.arcsin(np.cos(d(hoy)) * np.sin(np.radians(W(hoy))) / np.cos(ele(hoy)))

#%%
def thetai(hoy:np.array=HOYS_DEFAULT, inclination=90, azimuths=0): # incidence angle [in radians]

    g = np.degrees(azim(hoy)) - azimuths # if surface looks due S then azimuths=0
    return np.arccos(np.cos(ele(hoy)) * np.sin(np.radians(inclination)) * np.cos(np.radians(g)) 
        + np.sin(ele(hoy)) * np.cos(np.radians(inclination)))

def I0(hoy:np.array=HOYS_DEFAULT)->np.array: # Extra-terrestrial solar irradiance [W/m2]
    '''Appendix C in D.A. Katsaprakakis, Power Plant Synthesis, CRC Press, 2020.
    https://doi.org/10.1201/b22190.'''
    return 1373 * (1 + 0.033 * np.cos(np.radians(360*(hoy-3*24)/365*24)))

#%% ========================================== air mass
def air_mass(hoy:np.array=HOYS_DEFAULT, method:str= 'wiki' ):
    """wrapper function for the different air mass
    
    Args:
        hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
        method (str, optional): _description_. Defaults to 'wiki'.

    Returns:
        _type_: _description_
    """    
    
    dic = {'wiki': AM,
           'Kasten': AM2,
           'Kasten-Young': AM3,
           'Schoenberg' : AM4
           }
    if  method not in dic.keys():
        raise(ValueError(f"method can be [{dic.keys()}]"))
    return dic.get(method,None)(hoy)

def AM(hoy:np.array=HOYS_DEFAULT): # Air mass https://en.wikipedia.org/wiki/Air_mass_(solar_energy)
    AM =  1 / np.cos(z(hoy))
    return AM
def AM2(hoy:np.array=HOYS_DEFAULT):
    '''F. Kasten, A new table and approximation formula for the relative optical air mass, 
    Arch. Met. Geoph. Biokl. B. 14 (1965) 206–223. https://doi.org/10.1007/BF02248840.'''
    return 1 / (np.cos(z(hoy)) + 0.6556 * (6.379 - z(hoy))**-1.757)
def AM3(hoy:np.array=HOYS_DEFAULT):
    '''F. Kasten, A.T. Young, Revised optical air mass tables and approximation formula, 
    Appl. Opt., AO. 28 (1989) 4735–4738. https://doi.org/10.1364/AO.28.004735.'''
    return 1 / (np.cos(z(hoy)) + 0.50572 * (6.07995 - z(hoy))**-1.6364)
def AM4(hoy:np.array=HOYS_DEFAULT):
    '''E. Schoenberg, Theoretische Photometrie, in: K.F. Bottlinger, A. Brill, E. Schoenberg, 
    H. Rosenberg (Eds.), Grundlagen der Astrophysik, Springer, Berlin, Heidelberg, 1929: pp. 1–280. 
    https://doi.org/10.1007/978-3-642-90703-6_1.'''
    Re = 6371 # radius of the Earth [in km]
    yatm = 9 # effective height of the atmosphere [in km]
    r = Re / yatm
    return np.sqrt((r * np.cos(z(hoy)))**2 + 2 * r + 1) - r * np.cos(z(hoy))


#%% ======================== beam irradiance
def Ib(hoy:np.array=HOYS_DEFAULT): # Direct irradiance [in W/m2]
    '''https://www.pveducation.org/pvcdrom/properties-of-sunlight/air-mass
    The value of 1.353 kW/m2 is the solar constant and the number 0.7 arises from the fact that 
    about 70% of the radiation incident on the atmosphere is transmitted to the Earth. 
    The extra power term of 0.678 is an empirical fit to the observed data 
    and takes into account the non-uniformities in the atmospheric layers.
    A.B. Meinel, M.P. Meinel, Applied solar energy. An introduction, 
    Addison-Wesley Publishing Co.,Reading, MA, 1976. 
    https://www.osti.gov/biblio/7338398 (accessed January 7, 2021).'''
    a = np.power(abs(AM(hoy)), 0.678)
    return 1353 * np.power(0.7, a)

def Ib2(alt, hoy:np.array=HOYS_DEFAULT): 
    """Direct irradiance as function of height [in W/m2]

    '''E.G. Laue, The measurement of solar spectral irradiance at different terrestrial elevations, 
    Solar Energy. 13 (1970) 43–57. https://doi.org/10.1016/0038-092X(70)90006-X.'''

    Args:
        alt (_type_): altitude in ???? [m] or [km]
        hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.

    Returns:
        _type_: _description_
    """    
    a = np.power(np.abs(AM(hoy)), 0.678)
    return 1353 * ((1-0.14*alt) * np.power(0.7, a) + 0.14 * alt)

# def Ibs(hoy): #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<adapted by ASHRAE?
#     '''D.G. Stephenson, Equations for solar heat gain through windows, 
#     Solar Energy. 9 (1965) 81–86. https://doi.org/10.1016/0038-092X(65)90207-0.blank
#     C.A. Gueymard, Direct and indirect uncertainties in the prediction of tilted irradiance for solar engineering applications, 
#     Solar Energy. 83 (2009) 432–444. https://doi.org/10.1016/j.solener.2008.11.004.'''
#     return None

# 

# %%
