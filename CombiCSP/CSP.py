# -*- coding: utf-8 -*- 
"""
    @Author: G. Arnaoutakis, N. Papadakis
    @Date: 2022/  
"""

#TODO this module contains unused functions (not  tested)
#TODO    theta_transversal  function exists twice
#%%
'''Concentrating Solar Power plants'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

from iapws import IAPWS97
from CombiCSP import CtoK, HOYS_DEFAULT
import CombiCSP.SolarGeometry as sgh
from CombiCSP.SolarGeometry import W, z, d, thetai, azim, ele

#%%

#%%
# generator data

# SM = 2.5 # Solar Multiplier


#%% ===========================================================================
def theta_transversal(hoy:np.array=HOYS_DEFAULT)->float : 
    """Parabolic Trough theta  transversal incidence angle

    #TODO This function has the same name with another one in the same module
    
    #TODO  not tested
    

    Buscemi, A.; Panno, D.; Ciulla, G.; Beccali, M.; Lo Brano, V. 
    Concrete Thermal Energy Storage for Linear Fresnel Collectors: 
    Exploiting the South Mediterranean’s Solar Potential for Agri-Food Processes. 
    Energy Conversion and Management 2018, 166, 719–734, doi:10.1016/j.enconman.2018.04.075.

    Args:
        hoy (np.array): hour of year 

    Returns:
        float: theta  transversal incidence angle
    """    

    return np.arctan(np.sin(np.radians(azim(hoy))) * np.tan(np.radians(sgh.z(hoy))))


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
    return np.arctan(np.cos(np.radians(azim(hoy))) * np.tan(np.radians(sgh.z(hoy)))* np.cos(theta_transversal()))


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
    return np.arctan(abs(np.sin(azim(hoy)))/np.tan(ele(hoy)))

def thetai_longtitudinal(hoy:np.array=HOYS_DEFAULT): 
    return np.arcsin(np.cos(azim(hoy))*np.cos(ele(hoy)))


def shade_function(Ws,Wc, hoy:np.array=HOYS_DEFAULT):
    """Shade function for Parabolic Trough Solar Power Plants, 

    Stuetzle (2002) pp.29 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).
    N. Fraidenraich, C. Oliveira, A.F. Vieira da Cunha, J.M. Gordon, O.C. Vilela, 
    Analytical modeling of direct steam generation solar power plants, Solar Energy. 98 (2013) 511–522. 
    https://doi.org/10.1016/j.solener.2013.09.037.
    
    Args:
        Ws (_type_): #TODO  solar angle
        Wc (_type_): _description_
        hoy (_type_): _description_

    Returns:
        _type_: _description_
    """    
    return abs(Ws * np.cos(z(hoy)) / (Wc * np.cos(thetai(hoy))))

def end_loss(f,L,N, hoy:np.array=HOYS_DEFAULT):
    '''Lippke, 1995 in pp.31 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
    University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).'''
    return (1 - (f * np.tan(thetai(hoy)) / L)) * N

# def loss_regr(input_dict):
#     '''pp.36-42 in A.M. Patnode, Simulation and Performance Evaluation of Parabolic Trough Solar Power Plants, 
#     University of Wisconsin-Madison, 2006. https://minds.wisconsin.edu/handle/1793/7590 (accessed March 9, 2021).
#     '''
#     equation = a0 + a1 * T + a2 * T**2 + a3 * T**3 + Ib2(alt) * (b0 + b1 * T**2)
#     # Use dict flag to get {variable: value} output, not anonymous [value]
#     #solution = solve(equation.subs(input_dict), dict=True)
#     return equation

