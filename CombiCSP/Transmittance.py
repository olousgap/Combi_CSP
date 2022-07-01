# Currently this is an unused module
# Extensive testing needs to be carried out to validate the values.
#%%
import numpy as np

import CombiCSP.SolarGeometry as sgh
from CombiCSP.SolarGeometry import z

sin = np.sin
cos = np.cos
asin = np.arcsin
acos = np.arccos
atan = np.arctan
atan2 = np.arctan2
deg = np.degrees
rad = np.radians

HOYS_DEFAULT= sgh.HOYS_DEFAULT

#%%  ================================================= Transmittance
# the following functions relate to the [heliostats](https://www.nrel.gov/csp/heliocon.html) 
def Tr23km(alt, hoy:np.array=HOYS_DEFAULT): # Transmittance % Mid latitudes winter

    '''H.C. Hottel, A simple model for estimating the transmittance of direct solar radiation through clear atmospheres, 
    Solar Energy. 18 (1976) 129–134.'''
    a0 = 1.03 * (0.4327 - 0.00821 * (6 - alt)**2) #0.1283
    a1 = 1.01 * (0.5055 + 0.00595 * (6.5 - alt)**2) #0.7559
    k = 1.00 * (0.2711 + 0.01858 * (2.5 - alt)**2) #-0.3878
    return a0 + a1 * np.exp(-k)/np.exp(cos(np.radians(z(hoy)))) # needs rad despite z(hoy) already in rad???

def Tr5km(alt, hoy:np.array=HOYS_DEFAULT): # Transmittance % Mid latitudes winter
    '''H.C. Hottel, A simple model for estimating the transmittance of direct solar radiation through clear atmospheres, 
    Solar Energy. 18 (1976) 129–134.'''
    a0 = 1.04*(0.2538-0.0063*(6-alt)**2)
    a1 = 1.01*(0.7678+0.0010*(6.5-alt)**2)
    k = 1.00*(0.2490+0.0810*(2.5-alt)**2)
    return a0 + a1 * np.exp(-k)/np.exp(cos(np.radians(z(hoy)))) # needs rad despite z(hoy) already in rad???

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