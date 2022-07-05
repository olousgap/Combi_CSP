import pathlib
import numpy as np
import pandas as pd
import pytest
from CombiCSP import Economic_environment, SolarSystemLocation,  SolarTroughCalcs, HOYS_DEFAULT


@pytest.fixture
def stro()->SolarTroughCalcs:
    """SolarTrough ExampleData
    """    
    slobj = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)
    return SolarTroughCalcs(
        foc_len = 0.88 # [m] focal length CSPP T.1 in Mosleh19
        ,N = 1800 # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine
        ,L = 25 # [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine  
        ,Ws = 18 # [m] width between rows 18 INDITEP in pp.6 Fraidenraich13, pp.5 Zarza06
        ,Wr = 0.07 # tube outer diameter [m]
        ,Wc = 5.76,
        slobj= slobj
        )


@pytest.fixture
def ee():
    """Economic environment 
    """    
    return Economic_environment(   
            oil_price=60, 
            Eoil=11.766,
            currency_units='USD')
    
@pytest.fixture
def Ib():
    """SolarTower ExampleData
    """
    # pytest is configured in vscode at the root
    FNAME = pathlib.Path('tests/example_data/tmy_35.015_25.755_2005_2020.csv')
    pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
    Ib = pvgis.loc[:,'Gb(n)']
    return Ib
    
  

def test_ee_trough(stro, ee, Ib):
    o_tmp = stro.perform_calcs_NS(Ib, Tr=318, hoy=HOYS_DEFAULT)
    tmp_res_Dic = ee.economics_for_SolarTrough(oTr=o_tmp, 
        csp_area_costs =235, 
        csp_energy_price =248, 
        csp_discount_rate = 0.09, 
        power_block_cost = 910000.0,
        lifetime=range(30)
        )
    assert tmp_res_Dic['A_helio'] == 259200.0 
    np.testing.assert_almost_equal(
            tmp_res_Dic['cash_flow'][:2],
            [-117603680.04161312,   13795837.74])

    np.testing.assert_almost_equal(
        tmp_res_Dic['scenaria'][:8],
        (259200.0,
        82.28571428571428,
        62.29854949627815,
        74603.83,
        0.13670332645996,
        16.91753934467728,
        24129984.13280969,
        0.11252063804209933))
