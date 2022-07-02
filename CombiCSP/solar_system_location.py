# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
    @Credit: original functions based on G. Arnaoutakis
"""
import numpy as np
import pandas as pd
import pvlib

from CombiCSP import HOYS_DEFAULT
from CombiCSP.SolarGeometry import EoT,d 


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