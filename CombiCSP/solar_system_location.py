# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
    @Credit: original functions based on G. Arnaoutakis
"""
import numpy as np
import pandas as pd
import pvlib

HOYS_DEFAULT = np.arange(1, 8761, 1) # hours of year



class SolarSystemLocation:
    def __init__(self, lat:float,  lon:float, mer:float,dt_gmt:float, alt:float=0):
        """_summary_

        Args:
            lat (float): longitude of system 
            lon (float): longitude of system
            mer (float):  for Greece check to replace with 15 * dt_gmt (TODO better description)
            dt_gmt (float): time difference between Greenwich Mean Time
            alt (float): altitude of system (in m) (default: 0  - sea level)
        """        
        self.dt_gmt = dt_gmt # 
        self.lat = lat # Crete
        self.mer = mer # 
        self.lon = lon # Crete 35.2401° N, 24.8093° E [east negative, west positive

    @property
    def lat_rad(self)->float:
        """return latitude in radians

        Returns:
            _type_: _description_
        """        
        return np.radians(self.lat)
    @property	
    def long_rad(self)->float:
        """return longitude in radians

        Returns:
            float: _description_
        """        
        return np.radians(self.long)
    
    #%% ========================================== air mass
    def air_mass(self, hoy:np.array=HOYS_DEFAULT, method:str= 'wiki' ):
        """wrapper function for the different air mass
        
        Args:
            hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
            method (str, optional): _description_. Defaults to 'wiki'.

        Returns:
            _type_: _description_
        """    
        
        dic = {'wiki': self._AM_wiki,
            'Kasten': self._AM_Kasten,
            'Kasten-Young': self._AM3_KastenYoung,
            'Schoenberg' : self._AM_Shoenberg
            }
        if  method not in dic.keys():
            raise(ValueError(f"method can be [{dic.keys()}]"))
        return dic.get(method,None)(hoy)

    def _AM_wiki(self, hoy:np.array=HOYS_DEFAULT): # Air mass https://en.wikipedia.org/wiki/Air_mass_(solar_energy)
        AM =  1 / np.cos(self.z(hoy))
        return AM
    def _AM_Kasten(self, hoy:np.array=HOYS_DEFAULT):
        '''F. Kasten, A new table and approximation formula for the relative optical air mass, 
        Arch. Met. Geoph. Biokl. B. 14 (1965) 206–223. https://doi.org/10.1007/BF02248840.'''
        return 1 / (np.cos(self.z(hoy)) + 0.6556 * (6.379 - self.z(hoy))**-1.757)
    def _AM3_KastenYoung(self, hoy:np.array=HOYS_DEFAULT):
        '''F. Kasten, A.T. Young, Revised optical air mass tables and approximation formula, 
        Appl. Opt., AO. 28 (1989) 4735–4738. https://doi.org/10.1364/AO.28.004735.'''
        return 1 / (np.cos(self.z(hoy)) + 0.50572 * (6.07995 - self.z(hoy))**-1.6364)
    def _AM_Shoenberg(self, hoy:np.array=HOYS_DEFAULT):
        '''E. Schoenberg, Theoretische Photometrie, in: K.F. Bottlinger, A. Brill, E. Schoenberg, 
        H. Rosenberg (Eds.), Grundlagen der Astrophysik, Springer, Berlin, Heidelberg, 1929: pp. 1–280. 
        https://doi.org/10.1007/978-3-642-90703-6_1.'''
        Re = 6371 # radius of the Earth [in km]
        yatm = 9 # effective height of the atmosphere [in km]
        r = Re / yatm
        return np.sqrt((r * np.cos(self.z(hoy)))**2 + 2 * r + 1) - r * np.cos(self.z(hoy))

    def tsol(self, hoy:np.array=HOYS_DEFAULT): # solar time [in decimal hours] introduce if function for east/west<<<<<<<<<
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
        return hoy + (4*(self.lon-15*self.dt_gmt) + EoT(hoy))/60 # [60 min/h]


    def W(self, hoy:np.array=HOYS_DEFAULT): # solar hour angle [in degrees]
        #TODO rename to Omega
        '''given than for tsol = 12h it should be ω = 0ο and for the solar time range
        tsol = 0 – 24h the solar hour angle ranges from 0o to ±180'''
        return 15 * (self.tsol(hoy) - 12) # 360deg/24h = 15deg/h

    def ele(self, hoy:np.array=HOYS_DEFAULT): # solar elevation angle or solar height [in radians]
        return np.arcsin(np.cos(np.radians(self.lat)) * np.cos(d(hoy)) * np.cos(np.radians(self.W(hoy))) \
            + np.sin(np.radians(self.lat)) * np.sin(d(hoy)))

    def z(self, hoy:np.array=HOYS_DEFAULT):
        """Returns the solar zenith angle in radians

        solar zenith angle [in radians] https://en.wikipedia.org/wiki/Solar_zenith_angle   
        Args:
            hoy (np.array, optional): hour of year. Defaults to HOYS_DEFAULT.

        Returns:
        np.array: solar zenith angle in radians"""

        return np.arccos(np.cos(np.radians(self.lat)) * np.cos(d(hoy)) * np.cos(np.radians(self.W(hoy))) 
        + np.sin(np.radians(self.lat)) * np.sin(d(hoy)))

    def azim(self, hoy:np.array=HOYS_DEFAULT)->np.array: 
        """Returns the solar azimuth angle in radians

        Args:
            hoy (np.array, optional): hour of year. Defaults to HOYS_DEFAULT.

        Returns:
            np.array: solar azimuth angle in radians
        """   
        return np.arcsin(np.cos(d(hoy)) * np.sin(np.radians(self.W(hoy))) / np.cos(self.ele(hoy)))

# #TODO: consider creating a system/Unit parameters
#  def thetai(hoy:np.array=HOYS_DEFAULT, inclination=90, azimuths=0): # incidence angle [in radians]
#
#     g = deg(azim(hoy)) - azimuths # if surface looks due S then azimuths=0
#     return np.arccos(np.cos(ele(hoy)) * np.sin(np.radians(inclination)) * np.cos(np.radians(g)) 
#         + np.sin(ele(hoy)) * np.cos(np.radians(inclination)))

# ssCrete = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)


def get_pvgis_tmy_data(sysloc:SolarSystemLocation)->pd.DataFrame:
    """function that collects data from online. 

    #TODO this should be a class containing tmy data which can either be retrieved from a file on disk or PVGIS database
    

    Args:
        sysloc (_type_): _description_

    Returns:
        pd.DataFrame: _description_
    """    
    latitude = sysloc.lat
    longitude = sysloc.long
    OUTPUTFORMAT = 'json'

    dat = pvlib.iotools.get_pvgis_tmy(latitude, longitude, outputformat=OUTPUTFORMAT , usehorizon=True, userhorizon=None, 
        startyear=None, endyear=None, url='https://re.jrc.ec.europa.eu/api/', map_variables=True, timeout=30)

    df = dat[0]
    # dat[1] # monts of year for data 
    # dat[2] # metadata for the data set
    # dat[3] # variable explanation
    return df 

#%% Equations of time

def EoT(hoy:np.array=HOYS_DEFAULT): # equation of time [in minutes]
    gamma = 360*(hoy-1)/365
    return 2.2918*(0.0075+0.1868*np.cos(np.radians(gamma))-3.2077*np.sin(np.radians(gamma)) \
        -1.4615*np.cos(np.radians(2*gamma))-4.089*np.sin(np.radians(2*gamma)))

# TODO this do not work due to dayofyear --- uncomment when this is clear.
# def _calculate_simple_day_angle(dayofyear, offset=1):
#     """simple method for calculating the solar angle

#     Args:
#         dayofyear (_type_): _description_
#         offset (int, optional): _description_. Defaults to 1.

#     Returns:
#         _type_: solar angle  in radians 
#     """    
#     return (2. * np.pi / 365.) * (dayofyear - offset)

# def EoTS(hoy): # Equation of time from Duffie & Beckman and attributed to Spencer
#     #(1971) and Iqbal (1983) [in minutes]
#     day_angle = _calculate_simple_day_angle(dayofyear)
#     return (1440.0 / 2 / np.pi) * (0.0000075 +
#     0.001868 * np.cos(day_angle) - 0.032077 * np.sin(day_angle) -
#     0.014615 * np.cos(2.0 * day_angle) - 0.040849 * np.sin(2.0 * day_angle))

# def EoTPVCDROM(hoy): # equation of time [in minutes]
#     # PVCDROM: http://www.pveducation.org/pvcdrom/2-properties-sunlight/solar-time
#     # Soteris A. Kalogirou, "Solar Energy Engineering Processes and
#     # Systems, 2nd Edition" Elselvier/Academic Press (2009).
#     bday = _calculate_simple_day_angle(dayofyear) - (2.0 * np.pi / 365.0) * 80.0
#     return 9.87 * np.sin(2.0 * bday) - 7.53 * np.cos(bday) - 1.5 * np.sin(bday)


#%% ===================================== earth declination angles
def eda(hoy:np.array=HOYS_DEFAULT, method:str= 'wiki' ):
    """wrapper function for the different earth declination angles functions

    Useful reading [sunpos.py](https://levelup.gitconnected.com/python-sun-position-for-solar-energy-and-research-7a4ead801777)

    Args:
        hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.
        method (str, optional): _description_. Defaults to 'wiki'.

    Returns:
        _type_: _description_
    """    
    
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
    return 23.45 * np.sin(np.radians(360*(hoy+284*24)/365*24))

def d3(hoy:np.array=HOYS_DEFAULT):
    return 23.45 * np.sin(np.radians(360*(hoy-81*24)/365*24))

def d(hoy:np.array=HOYS_DEFAULT): # [in radians] https://en.wikipedia.org/wiki/Sunrise_equation
    '''Calculate ecliptic coordinates (ecliptic longitude and obliquity of the
    ecliptic in radians but without limiting the angle to be less than 2*Pi
    (i.e., the result may be greater than 2*Pi)'''
    dOmega = 2.1429 - 0.0010394594 * hoy
    dMeanLongitude = 4.8950630 + 0.017202791698 * hoy
    dMeanAnomaly = 6.2400600 + 0.0172019699 * hoy
    dEclipticLongitude = dMeanLongitude + 0.03341607 * np.sin(dMeanAnomaly) 
    + 0.00034894 * np.sin( 2 * dMeanAnomaly) - 0.0001134 - 0.0000203 * np.sin(dOmega)
    dEclipticObliquity = 0.4090928 - 6.2140e-9 * hoy 
    + 0.0000396 * np.cos(dOmega)
    '''Calculate celestial coordinates ( right ascension and declination ) in radians
    but without limiting the angle to be less than 2*Pi (i.e., the result may be
    greater than 2*Pi)'''
    dSin_EclipticLongitude = np.sin( dEclipticLongitude )
    dY = np.cos( dEclipticObliquity ) * dSin_EclipticLongitude
    dX = np.cos( dEclipticLongitude )
    dRightAscension = np.arctan2( dY,dX )
    if dRightAscension.any() < 0: dRightAscension = dRightAscension + 2*np.pi
    return np.arcsin( np.sin( dEclipticObliquity ) * dSin_EclipticLongitude )

