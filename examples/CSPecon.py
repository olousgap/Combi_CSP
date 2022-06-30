#%%
from pylab import *
import numpy as np
import numpy_financial as npf
import pandas as pd 
from CombiCSP.SolarGeometry import *
import CombiCSP.CSP as CSP
from CSPCret import *

'''Engineering inputs from CSPCret'''
Pcsp = Ptrough # 250 [kW] Mosleh19
Ecsp = Etrough # 1713.200 [annual MWh] Mosleh19

#Ac = 3720 # replace with Ac(Wc, L, N) from CSPCret

mpcm = pcm_rho * 200 #mpcm = pcm_rho * Vpcm # kg/m3 * m3

#Paux = 9 # [MW] T.5 Pantaleo17
Eaux = 20 # [MWh]
Eoil = Eaux*0.5883 # [BOE] 1MWh = 0.5883BOE https://www.convert-me.com/en/convert/energy/kwh/kwh-to-boe.html?u=kwh&v=1%2C000
Egas = Eaux*3.412 # [m BTU] 1MWh = 3.412mBTU https://www.convert-me.com/en/convert/energy/kwh/kwh-to-mymmbtu.html?u=kwh&v=1%2C000

def mbtu_m3(mbtu:float):
    """Converts MMBtu (1 Million BTU units) to m^3 of natural gas

    Args:
        mbtu (_type_): _description_

    Returns:
        _type_: _description_
    """    
    '''
    https://www.indexmundi.com/commodities/glossary/mmbtu#:~:text=Natural%20gas%20is%20measured%20in,on%20quality%2C%20when%20burned).
    '''
    return mbtu * 28.263682

'''Financial inputs'''
#discount_rate = 0.1 # Mosleh19 0.06...0.1

'''[energy_price_E_MWh],[discount_rate] Υπουργική Απόφαση ΥΠΕΝ/ΔΑΠΕΕΚ/30971/1190/2020 - ΦΕΚ 1045/Β/26-3-2020'''
energy_price = [[248,268,176,153,185,133,90,87,87],
                [0.09,0.09,0.08,0.074,0.08,0.074,0.08,0.08,0.074]]
energy_price = pd.DataFrame(energy_price,
    columns = ['CSP','CSP+Storage2h','Biomass 1MW','Biomass 5MW',
    'Biogas 1MW','Biogas 5MW','Hydro 3MW','Hydro 15MW','PV6KWroof'])

# energy_price = pd.DataFrame( [[248,268,176,153,185,133,90,87,87],
#                 [0.09,0.09,0.08,0.074,0.08,0.074,0.08,0.08,0.074]],
#     columns = ['CSP','CSP+Storage2h','Biomass 1MW','Biomass 5MW',
#     'Biogas 1MW','Biogas 5MW','Hydro 3MW','Hydro 15MW','PV6KWroof'], index= ['price', 'discount_rate']).transpose()



csp_energy_price = energy_price['CSP'].loc[0]
csp_pcm_energy_price = energy_price['CSP+Storage2h'].loc[0]
csp_discount_rate = energy_price['CSP'].loc[1]
bio_energy_price = energy_price['Biomass 5MW'].loc[0]
biogas_energy_price = energy_price['Biogas 5MW'].loc[0]
bio_discount_rate = energy_price['Biomass 5MW'].loc[1]

oil_price = 60 # np.arange(12, 112, 10)# [$/barrel] https://www.statista.com/statistics/262860/uk-brent-crude-oil-price-changes-since-1976/
gas_price = 7 # np.arange(2.5, 7.1, 2)# [$/m BTU/year] https://www.statista.com/statistics/252791/natural-gas-prices/

# alternative dataframes
#data = [[248,0.09],[268,0.09],[176,0.08],[153,0.074],[185,0.08],[133,0.074],[90,0.08],[87,0.08],[87,0.074]]
#[['csp',248,0.09],['csp_stor2h',268,0.09],['bio_burn_1MW',176,0.08],['bio_burn_5MW',153,0.074],['bio_gas1MW',185,0.08],['bio_gas5MW',133,0.074],['hydro3MW',90,0.08],['hydro15MW',87,0.08],['pv6KWroof',87,0.074]]
#df = pd.DataFrame(energy_price, columns = ['Name', 'price (E/MWh)','discount_rate (%)']).transpose()
#df[1].iloc[1]

'''T.4 Turchi19 all in units [$/MWh]'''
energy_cost = [6.2e4]
#energy_cost = pd.DataFrame(power_cost, columns = ['thermal_energy_storage'])

'''all in units [$/kg] U. Herrmann, D.W. Kearney, 
Survey of Thermal Energy Storage for Parabolic Trough Power Plants, 
Journal of Solar Energy Engineering. 124 (2002) 145–152. 
https://doi.org/10.1115/1.1467601.'''
# taken from CSP.py pcm dataframe
# pcm_cost = pcm['<salt>'].loc[6]

'''all in units [$/MW] T.4 Turchi19, stokerCHP in pp.17 Biomass for Heat and Power IRENA report, 1e6/9.050=1.1e5 euro/MWth T.7 Pantaleo17
# T.4 Turchi19, 1.16e6 T.2 Turchi10, 2.2e6 T.7 Pantaleo17'''
power_cost = [[9.1e5, 0.9e5, 4e6, 6e6]]
power_cost = pd.DataFrame(
    power_cost,
    columns = ['power_block_cost','balance_of_plant','boiler_cost','gasifier_cost'])
boiler_cost = power_cost['boiler_cost'].loc[0]
gasifier_cost = power_cost['gasifier_cost'].loc[0]
power_block_cost = power_cost['power_block_cost'].loc[0]

'''T.4 Turchi19 all in units [$/m2]'''
csp_area_costs = np.sum([25, 150, 60])
csp_area_costs_df = [[25, 150, 60]]
csp_area_costs_df = pd.DataFrame(
    csp_area_costs_df, index=['USD/m${^2}$'],
    columns = ['site_dev','coll_cost','htf_cost'])

lifetime = np.arange(0, 31, 1)

# %%
