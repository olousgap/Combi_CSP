#%%
'''Geometry calculations as a function of hourly local time     ALT + SHIFT +0 to unfold levels
                                                                ALT + 0 to fold levels'''
import numpy as np
# from pylab import *
#for ineichen
import pandas as pd
import pvlib
from pvlib import clearsky, atmosphere, solarposition

sin = np.sin
cos = np.cos
asin = np.arcsin
acos = np.arccos
atan = np.arctan
atan2 = np.arctan2
deg = np.degrees
rad = np.radians

# args = hoy = np.arange(1, 8761, 1) # hour of year
HOYS_DEFAULT = np.arange(1, 8761, 1) # hours of year


#TODO Move this parameters to another location because now the system can only calculate for these long and lat.
dt_gmt = +2 # time difference between Greenwich Mean Time
lat = 35 # Crete
mer = -25 # for Greece check to replace with 15 * dt_gmt
lon = 24 # Crete 35.2401° N, 24.8093° E [east negative, west positive]

#%% ===================================== earth declination angles
def eda(hoy:np.array=HOYS_DEFAULT, method= 'wiki' ):
    dic = {'wiki':d,
           'Katsaprakakis': d2,
           '-81': d3,
           'pveducation' :d1
           }
    if  method not in dic.keys():
        raise(ValueError(f"method can be [{dic.keys()}]"))
    return dic.get(method,None)(hoy)

def d1(hoy:np.array=HOYS_DEFAULT): # earth declination angle [in degrees]
    """earth declination angle 

    The +10 comes from the fact that the winter solstice occurs before the start
    of the year. The equation also assumes that the suns orbit is a perfect circle and 
    the factor of 360/365 converts the day number to a position in the orbit.

    Args:
        hoy (np.array, optional): an array in hours of year . Defaults to HOYS_DEFAULT.

    Returns:
        _type_: earth inclination angle in degrees
    """
    # https://www.pveducation.org/pvcdrom/properties-of-sunlight/declination-angle#footnote1_osno74c'''
    return -23.45 * np.cos( np.radians(360*(hoy+10*24)/365*24))

def d2(hoy:np.array=HOYS_DEFAULT):
    '''pp.6, Appendix C in D.A. Katsaprakakis, Power Plant Synthesis, CRC Press, 2020.
    https://doi.org/10.1201/b22190.'''
    return 23.45 * sin(rad(360*(hoy+284*24)/365*24))

def d3(hoy:np.array=HOYS_DEFAULT):
    return 23.45 * sin(rad(360*(hoy-81*24)/365*24))

def d(hoy:np.array=HOYS_DEFAULT): # [in radians] https://en.wikipedia.org/wiki/Sunrise_equation
    '''Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
    ecliptic in radians but without limiting the angle to be less than 2*Pi
    (i.e., the result may be greater than 2*Pi)'''
    dOmega = 2.1429 - 0.0010394594 * hoy
    dMeanLongitude = 4.8950630 + 0.017202791698 * hoy
    dMeanAnomaly = 6.2400600 + 0.0172019699 * hoy
    dEclipticLongitude = dMeanLongitude + 0.03341607 * sin(dMeanAnomaly) 
    + 0.00034894 * sin( 2 * dMeanAnomaly) - 0.0001134 - 0.0000203 * sin(dOmega)
    dEclipticObliquity = 0.4090928 - 6.2140e-9 * hoy 
    + 0.0000396 * cos(dOmega)
    '''Calculate celestial coordinates ( right ascension and declination ) in radians
    but without limiting the angle to be less than 2*Pi (i.e., the result may be
    greater than 2*Pi)'''
    dSin_EclipticLongitude = sin( dEclipticLongitude )
    dY = cos( dEclipticObliquity ) * dSin_EclipticLongitude
    dX = cos( dEclipticLongitude )
    dRightAscension = atan2( dY,dX )
    if dRightAscension.any() < 0: dRightAscension = dRightAscension + 2*np.pi
    return asin( sin( dEclipticObliquity ) * dSin_EclipticLongitude )

#%% Equations of time

def EoT(hoy:np.array=HOYS_DEFAULT): # equation of time [in minutes]
    gamma = 360*(hoy-1)/365
    return 2.2918*(0.0075+0.1868*cos(rad(gamma))-3.2077*sin(rad(gamma)) \
        -1.4615*cos(rad(2*gamma))-4.089*sin(rad(2*gamma)))

def _calculate_simple_day_angle(dayofyear, offset=1):
    """simple method for calculating the solar angle

    Args:
        dayofyear (_type_): _description_
        offset (int, optional): _description_. Defaults to 1.

    Returns:
        _type_: solar angle  in radians 
    """    
    return (2. * np.pi / 365.) * (dayofyear - offset)

def EoTS(hoy): # Equation of time from Duffie & Beckman and attributed to Spencer
    #(1971) and Iqbal (1983) [in minutes]
    day_angle = _calculate_simple_day_angle(dayofyear)
    return (1440.0 / 2 / np.pi) * (0.0000075 +
    0.001868 * np.cos(day_angle) - 0.032077 * np.sin(day_angle) -
    0.014615 * np.cos(2.0 * day_angle) - 0.040849 * np.sin(2.0 * day_angle))

def EoTPVCDROM(hoy): # equation of time [in minutes]
    # PVCDROM: http://www.pveducation.org/pvcdrom/2-properties-sunlight/solar-time
    # Soteris A. Kalogirou, "Solar Energy Engineering Processes and
    # Systems, 2nd Edition" Elselvier/Academic Press (2009).
    bday = _calculate_simple_day_angle(dayofyear) - (2.0 * np.pi / 365.0) * 80.0
    return 9.87 * np.sin(2.0 * bday) - 7.53 * np.cos(bday) - 1.5 * np.sin(bday)

#%% 
def tsol(hoy:np.array=HOYS_DEFAULT): # solar time [in decimal hours] introduce if function for east/west<<<<<<<<<
    '''CALCULATION NEGLECTS DAYLIGHT SAVING ON SUMMER
    https://www.pveducation.org/pvcdrom/properties-of-sunlight/the-suns-position
    hoy + ((lat-mer)/15 + EoT(hoy))/60 # [60 min/h]
    pp.5, Appendix C in D.A. Katsaprakakis, Power Plant Synthesis, CRC Press, 2020.
    https://doi.org/10.1201/b22190'''
    return hoy + (4*(lon-15*dt_gmt) + EoT(hoy))/60 # [60 min/h]


def W(hoy:np.array=HOYS_DEFAULT): # solar hour angle [in degrees]
    '''given than for tsol = 12h it should be ω = 0ο and for the solar time range
    tsol = 0 – 24h the solar hour angle ranges from 0o to ±180'''
    return 15 * (tsol(hoy) - 12) # 360deg/24h = 15deg/h

def ele(hoy:np.array=HOYS_DEFAULT): # solar elevation angle or solar height [in radians]
    return asin(cos(rad(lat)) * cos(d(hoy)) * cos(rad(W(hoy)))
    + sin(rad(lat)) * sin(d(hoy)))

def z(hoy:np.array=HOYS_DEFAULT): # solar zenith angle [in radians] https://en.wikipedia.org/wiki/Solar_zenith_angle
    return acos(cos(rad(lat)) * cos(d(hoy)) * cos(rad(W(hoy))) 
    + sin(rad(lat)) * sin(d(hoy)))

def azim(hoy:np.array=HOYS_DEFAULT): # solar azimuth [in radians]
    return asin(cos(d(hoy)) * sin(rad(W(hoy))) / cos(ele(hoy)))

#%%
def thetai(hoy:np.array=HOYS_DEFAULT): # incidence angle [in radians]
    inclination=90
    azimuths=0
    g = deg(azim(hoy)) - azimuths # if surface looks due S then azimuths=0
    return acos(cos(ele(hoy)) * sin(rad(inclination)) * cos(rad(g)) 
    + sin(ele(hoy)) * cos(rad(inclination)))

def I0(hoy:np.array=HOYS_DEFAULT): # Extra-terrestrial solar irradiance [W/m2]
    '''Appendix C in D.A. Katsaprakakis, Power Plant Synthesis, CRC Press, 2020.
    https://doi.org/10.1201/b22190.'''
    return 1373 * (1 + 0.033 * cos(rad(360*(hoy-3*24)/365*24)))

#%% ========================================== air mass
def AM(hoy:np.array=HOYS_DEFAULT): # Air mass https://en.wikipedia.org/wiki/Air_mass_(solar_energy)
    AM =  1 / cos(z(hoy))
    return AM
def AM2(hoy:np.array=HOYS_DEFAULT):
    '''F. Kasten, A new table and approximation formula for the relative optical air mass, 
    Arch. Met. Geoph. Biokl. B. 14 (1965) 206–223. https://doi.org/10.1007/BF02248840.'''
    return 1 / (cos(z(hoy)) + 0.6556 * (6.379 - z(hoy))**-1.757)
def AM3(hoy:np.array=HOYS_DEFAULT):
    '''F. Kasten, A.T. Young, Revised optical air mass tables and approximation formula, 
    Appl. Opt., AO. 28 (1989) 4735–4738. https://doi.org/10.1364/AO.28.004735.'''
    return 1 / (cos(z(hoy)) + 0.50572 * (6.07995 - z(hoy))**-1.6364)
def AM4(hoy:np.array=HOYS_DEFAULT):
    '''E. Schoenberg, Theoretische Photometrie, in: K.F. Bottlinger, A. Brill, E. Schoenberg, 
    H. Rosenberg (Eds.), Grundlagen der Astrophysik, Springer, Berlin, Heidelberg, 1929: pp. 1–280. 
    https://doi.org/10.1007/978-3-642-90703-6_1.'''
    Re = 6371 # radius of the Earth [in km]
    yatm = 9 # effective height of the atmosphere [in km]
    r = Re / yatm
    return np.sqrt((r * cos(z(hoy)))**2 + 2 * r + 1) - r * cos(z(hoy))


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
    a = np.power(abs(AM(hoy)), 0.678)
    return 1353 * ((1-0.14*alt) * np.power(0.7, a) + 0.14 * alt)
# def Ibs(hoy): #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<adapted by ASHRAE?
#     '''D.G. Stephenson, Equations for solar heat gain through windows, 
#     Solar Energy. 9 (1965) 81–86. https://doi.org/10.1016/0038-092X(65)90207-0.blank
#     C.A. Gueymard, Direct and indirect uncertainties in the prediction of tilted irradiance for solar engineering applications, 
#     Solar Energy. 83 (2009) 432–444. https://doi.org/10.1016/j.solener.2008.11.004.'''
#     return None

#%%  ================================================= Transmittance
# the following functions relate to the [heliostats](https://www.nrel.gov/csp/heliocon.html) 
def Tr23km(alt, hoy:np.array=HOYS_DEFAULT): # Transmittance % Mid latitudes winter

    '''H.C. Hottel, A simple model for estimating the transmittance of direct solar radiation through clear atmospheres, 
    Solar Energy. 18 (1976) 129–134.'''
    a0 = 1.03 * (0.4327 - 0.00821 * (6 - alt)**2) #0.1283
    a1 = 1.01 * (0.5055 + 0.00595 * (6.5 - alt)**2) #0.7559
    k = 1.00 * (0.2711 + 0.01858 * (2.5 - alt)**2) #-0.3878
    return a0 + a1 * np.exp(-k)/np.exp(cos(rad(z(hoy)))) # needs rad despite z(hoy) already in rad???

def Tr5km(alt, hoy:np.array=HOYS_DEFAULT): # Transmittance % Mid latitudes winter
    '''H.C. Hottel, A simple model for estimating the transmittance of direct solar radiation through clear atmospheres, 
    Solar Energy. 18 (1976) 129–134.'''
    a0 = 1.04*(0.2538-0.0063*(6-alt)**2)
    a1 = 1.01*(0.7678+0.0010*(6.5-alt)**2)
    k = 1.00*(0.2490+0.0810*(2.5-alt)**2)
    return a0 + a1 * np.exp(-k)/np.exp(cos(rad(z(hoy)))) # needs rad despite z(hoy) already in rad???

def TrD23km(R): 
    """transmmitance at 23km 

    Args:
        R (_type_): _description_

    Returns:
        _type_: _description_
    """    
    return 0.6739+10.46*R-1.7*R**2+0.2845*R**3

def TrD5km(R): 
    """transmmitance at 23km 

    Args:
        R (_type_): _description_

    Returns:
        _type_: _description_
    """    
    return 1.293+27.48*R-3.394*R**2+0*R**3

def TrV23km(R): 
    return 0.99326-0.1046*R+0.017*R**2-0.002845*R**3

def TrV5km(R): 
    """_summary_

    Args:
        R (_type_):  slant range in [km]

    Returns:
        _type_: _description_
    """    
    return 0.98707-0.2748*R+0.03394*R**2-0*R**3


def TrVH(Ht,R,alt): #V_H transmittance model [all units in km]
    '''C.L. Pitman, L.L. Vant-Hull, Atmospheric transmittance model for a solar beam propagating between a heliostat and a receiver, 
    Sandia National Labs., Albuquerque, NM (USA); Houston Univ., TX (USA). Energy Lab., 1984. https://doi.org/10.2172/5148541.'''
    #Ht = 0.100 # Tower height [km]
    b = 0.17 # attenuation coefficient due to scattering by aerosols and air molecules at 550nm [km-1]
    rw = 5.9 # water vapor density [g/m3]
    '''G.P. Anderson, et al , Reviewing atmospheric radiative transfer modeling: 
    new developments in high- and moderate-resolution FASCODE/FASE and MODTRAN, in: 
    Optical Spectroscopic Techniques and Instrumentation for Atmospheric and Space Research II, 
    International Society for Optics and Photonics, 1996: pp. 82–93. https://doi.org/10.1117/12.256105.'''
    A = (0.0112 * alt + 0.0822) * np.log(((b + 0.0003 * rw)/0.00455)) # tower focal height [km-1]
    S = 1 - ((0.00101 * rw + 0.0507) * np.sqrt(b + 0.0091))
    # necessary because the solar beam is composed of many wavelengths and because absorption 
    #(which is not described by b) also occurs along the path. (In the case of single-wavelength transmission, S is always unity.)
    C = (0.0105 * rw + 0.724) * (b - 0.0037)**S
    ks = C * np.exp(-A * Ht) # broadband extinction coefficient between tower and receiver atop [km-S]
    return np.exp(-ks * R**S) # R slant range in [km]

# %%
