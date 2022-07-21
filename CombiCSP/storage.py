# -*- coding: utf-8 -*- 
"""
    @Author: G. Arnaoutakis, N. Papadakis
    @Date: 2022/mm/dd
    @Credit: original functions from G. Arnaoutakis
"""

""" 
    This module contains 
    - data relating to phase change materials 
    - functions relating to heat energy storage

    TODO: This module should be reviewed and rewritten.
"""


#%%
'''Concentrating Solar Power plants'''
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from scipy import integrate

import CombiCSP.SolarGeometry as sgh
# from CombiCSP.SolarGeometry import z,d, W
# from CombiCSP.SolarGeometry import thetai, azim, ele
# from iapws import IAPWS97
# from numpy import pi, sin, cos, tan
#%% material properties for storage
#PCM data  Phase Change Material
pcm_data = [[2200,2257,2110,2044,2380],
            [212,174,226,149.7,280],
            [282,308,333,380,250,802],
            [np.NaN,0.5,0.5,0.5,0.52,5.0],
            [1.733,1.588,1.240,np.NaN,np.NaN,np.NaN],
            [2.553,1.650,1.341,np.NaN,np.NaN,np.NaN],
            [0.2,0.2,0.3,1.0,np.NaN,0.15]] #NaNO2 assumption
pcm = pd.DataFrame(pcm_data,
    columns = ['NaNO2','NaNO3','KNO3','KOH','H250','NaCl']
    #, index = ['density', 'latent_heat', 'melting_point', 'thermal conduct. coeff', 'c_p Solid', 'c_p Liquid', 'cost' ]
    )
pcm_rho = pcm['NaNO3'].loc[0]
pcm_latent_heat = pcm['NaNO3'].loc[1]#*0.27778 # [kJ/kgK to Wh/kgK] https://www.cactus2000.de/uk/unit/masscp1.php
pcm_melting_point = pcm['NaNO3'].loc[2] # [oC]
pcm_thermal_conductivity_coeff = pcm['NaNO3'].loc[3]
pcm_solid_cp = pcm['NaNO3'].loc[4] # [kJ/kgK]
pcm_liquid_cp = pcm['NaNO3'].loc[5] # [kJ/kgK]
pcm_cost = pcm['NaNO3'].loc[6]


# define materials
Vpcm = 200 # [m3] Volume
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
# T = CSCP(Tr)



#%% storage

def heatloss_Mertins(d,epsilon,dt):
    """ the heat loss from the receiver

    Args:
        d (_type_): diameter (m)
        epsilon (_type_): emmissivity  []
        dt (_type_): the temperature difference between the tube and ambient [K]

    Returns:
        _type_: _description_
    """    
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
    return R * np.log(x_a)

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


# def pcm_temp():
#     #TODO this function did not work (T, cpcm were not available)
#     '''https://stackabuse.com/solving-systems-of-linear-equations-with-pythons-numpy/
#     Y.-Q. Li, Y.-L. He, Z.-F. Wang, C. Xu, W. Wang, Exergy analysis of two phase change materials storage system 
#     for solar thermal power with finite-time thermodynamics, Renewable Energy. 39 (2012) 447–454. 
#     https://doi.org/10.1016/j.renene.2011.08.026.
#     H.J. Mosleh, R. Ahmadi, Linear parabolic trough solar power plant assisted with latent thermal energy storage system: 
#     A dynamic simulation, Applied Thermal Engineering. 161 (2019) 114204. https://doi.org/10.1016/j.applthermaleng.2019.114204.
#     '''
#     Tin = Tr
#     Nc = hp * Ap / mpcm * cpcm
#     T2 = Tmelt + (T - Tmelt) * np.exp(-Nc)
#     return T2 #Estore

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
