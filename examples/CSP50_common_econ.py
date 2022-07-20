# This is a file that declares common economic parameters and data for CSP50 analysis. 
# 
# It is created separately during developemt, with the idea being to mimisize the code and avoid duplication between the different versions of CSP50
#  All the CSP50 variations should import this file - until the code base is mature enough.
#%%
import numpy as np
import pandas as pd
from CombiCSP.storage import pcm_rho

#%%

'''Engineering inputs from CSPCret'''
# Pcsp = Ptrough # 250 [kW] Mosleh19
# Ecsp = Etrough # 1713.200 [annual MWh] Mosleh19
# Ac = 3720 # replace with Ac(Wc, L, N) from CSPCret

Vpcm = 200            # Volume Phase change materials
mpcm = pcm_rho *Vpcm  # mass of phase change material
                      # mpcm = pcm_rho * Vpcm # kg/m3 * m3

#Paux = 9 # [MW] T.5 Pantaleo17
Eaux = 20 # [MWh]
Eoil = Eaux*0.5883 # [BOE] 1MWh = 0.5883BOE https://www.convert-me.com/en/convert/energy/kwh/kwh-to-boe.html?u=kwh&v=1%2C000
Egas = Eaux*3.412 # [m BTU] 1MWh = 3.412mBTU https://www.convert-me.com/en/convert/energy/kwh/kwh-to-mymmbtu.html?u=kwh&v=1%2C000

'''Financial inputs'''
#discount_rate = 0.1 # Mosleh19 0.06...0.1


''' The followind data pertain to Greece based on:
[energy_price_E_MWh],[discount_rate] Υπουργική Απόφαση ΥΠΕΝ/ΔΑΠΕΕΚ/30971/1190/2020 - ΦΕΚ 1045/Β/26-3-2020'''
energy_price_GR_df = pd.DataFrame( [
    ['CSP' , 248, 0.09],
    ['CSP+Storage2h', 268,  0.09],
    ['Biomass 1MW',   176,  0.08],
    ['Biomass 5MW',   153,  0.074],
    ['Biogas 1MW',    185,  0.08],
    ['Biogas 5MW',    133,  0.074],
    ['Hydro 3MW',      90,  0.08],
    ['Hydro 15MW',     87,  0.08],
    ['PV6KWroof',      87,  0.074]
    ],
    columns = ['type', 'price', 'discount_rate']).set_index('type', drop=True)
print(energy_price_GR_df)

csp_energy_price = energy_price_GR_df.loc['CSP','price']
csp_pcm_energy_price = energy_price_GR_df.loc['CSP+Storage2h', 'price']
csp_discount_rate = energy_price_GR_df.loc['CSP','discount_rate']

bio_energy_price = energy_price_GR_df.loc['Biomass 5MW', 'price']
biogas_energy_price = energy_price_GR_df.loc['Biogas 5MW', 'price']
bio_discount_rate = energy_price_GR_df.loc['Biomass 5MW','discount_rate']

oil_price = 60 # np.arange(12, 112, 10)# [$/barrel] https://www.statista.com/statistics/262860/uk-brent-crude-oil-price-changes-since-1976/
gas_price = 7 # np.arange(2.5, 7.1, 2)# [$/m BTU/year] https://www.statista.com/statistics/252791/natural-gas-prices/

# alternative dataframes
#data = [[248,0.09],[268,0.09],[176,0.08],[153,0.074],[185,0.08],[133,0.074],[90,0.08],[87,0.08],[87,0.074]]
#[['csp',248,0.09],['csp_stor2h',268,0.09],['bio_burn_1MW',176,0.08],['bio_burn_5MW',153,0.074],['bio_gas1MW',185,0.08],['bio_gas5MW',133,0.074],['hydro3MW',90,0.08],['hydro15MW',87,0.08],['pv6KWroof',87,0.074]]
#df = pd.DataFrame(energy_price, columns = ['Name', 'price (E/MWh)','discount_rate (%)']).transpose()
#df[1].iloc[1]



'''T.4 Turchi19 all in units [$/MWh]'''
energy_cost = {'thermal_energy_storage': 6.2e4}

'''all in units [$/kg] U. Herrmann, D.W. Kearney, 
Survey of Thermal Energy Storage for Parabolic Trough Power Plants, 
Journal of Solar Energy Engineering. 124 (2002) 145–152. 
https://doi.org/10.1115/1.1467601.'''
# taken from CSP.py pcm dataframe
# pcm_cost = pcm['<salt>'].loc[6]

'''all in units [$/MW] T.4 Turchi19, stokerCHP in pp.17 Biomass for Heat and Power IRENA report, 1e6/9.050=1.1e5 euro/MWth T.7 Pantaleo17
# T.4 Turchi19, 1.16e6 T.2 Turchi10, 2.2e6 T.7 Pantaleo17'''

dic_Power_cost = {
    'power_block_cost': 9.1e5,
    'balance_of_plant':0.9e5,
    'boiler_cost': 4e6,
    'gasifier_cost': 6e6
}

boiler_cost = dic_Power_cost['boiler_cost']
gasifier_cost = dic_Power_cost['gasifier_cost']
power_block_cost = dic_Power_cost['power_block_cost']

'''T.4 Turchi19 all in units [$/m2]'''

csp_area_costs_df = pd.DataFrame(
    [[25, 150, 60]], index=['USD/m${^2}$'],
    columns = ['site_dev','coll_cost','htf_cost'])
csp_area_costs = csp_area_costs_df.iloc[0,:].sum()

lifetime = np.arange(0, 31, 1)
