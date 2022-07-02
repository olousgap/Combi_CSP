#TODO this module contains unused functions
#     they have not been tested and in some cases there are issues
#     eg. theta_transversal exists twice
#%%
'''Concentrating Solar Power plants'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

from iapws import IAPWS97
from CombiCSP import CtoK
import CombiCSP.SolarGeometry as sgh
from CombiCSP.SolarGeometry import W, z, d, thetai, azim, ele, HOYS_DEFAULT

#%%
# from numpy import log as ln
# NaN = np.nan
# pi = np.pi
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



#%%
# generator data

# SM = 2.5 # Solar Multiplier


#%% ===========================================================================

def theta_i(hoy:np.array=HOYS_DEFAULT)->float: 
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
    return arctan(cos(rad(azim(hoy))) * tan(rad(sgh.z(hoy)))* cos(theta_transversal()))



def theta_transversal(hoy:np.array=HOYS_DEFAULT)->float : 
    """Parabolic Trough theta  transversal incidence angle

    #TODO This function has the same name with another one in the same module
    

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

    return np.arctan(sin(rad(azim(hoy))) * tan(rad(sgh.z(hoy))))



# not tested
def thetai_transversal(hoy:np.array=HOYS_DEFAULT):
    """_summary_

    #TODO This function has the same name with another one in the same module
    
    Morin, G.; Dersch, J.; Platzer, W.; Eck, M.; Häberle, A. 
    Comparison of Linear Fresnel and Parabolic Trough Collector Power Plants. 
    Solar Energy 2012, 86, 1–12, doi:10.1016/j.solener.2011.06.020.

    
    Args:
        hoy (np.array, optional): _description_. Defaults to HOYS_DEFAULT.

    Returns:
        _type_: _description_
    """    
    return np.arctan(abs(sin(azim(hoy)))/tan(ele(hoy)))

def thetai_longtitudinal(hoy:np.array=HOYS_DEFAULT): 
    return np.arcsin(cos(azim(hoy))*cos(ele(hoy)))


def shade_function(Ws,Wc, hoy:np.array=HOYS_DEFAULT):
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
    return abs(Ws * cos(sgh.z(hoy)) / (Wc * cos(thetai(hoy))))

def end_loss(f,L,N, hoy:np.array=HOYS_DEFAULT):
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

