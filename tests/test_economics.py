import pathlib
import numpy as np
import pandas as pd
import pytest
from CombiCSP import Economic_environment, SolarSystemLocation, SolarTowerCalcs


@pytest.fixture
def st():
    """SolarTower ExampleData
    """    
    slobj = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)
    return SolarTowerCalcs(alt = 200*10e-3 , Ht = 0.1, Ar = 99.3 , A_helio = 225000, slobj=slobj)

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
    

def test_ee(ee):
    assert ee.oil_price == 60
    assert ee._Eoil == 11.766
    
def test_ee_tower(st, ee, Ib):
    o_tmp = st.perform_calc(Ib)
    tmp_res_Dic =ee.economics_for_Solar_tower(
            oTow= o_tmp,
            csp_area_costs= 235,
            csp_energy_price=248,
            csp_discount_rate= 0.09,
            power_block_cost=910000.0,
            capital_csp=5000000,
        lifetime=range(30))
    assert tmp_res_Dic['A_helio'] == 225000
    np.testing.assert_almost_equal(
            tmp_res_Dic['cash_flow_tow'][:2],
            [-106352083.86326265, 23230996.15] )

    np.testing.assert_almost_equal(
        tmp_res_Dic['tow_scenaria'][:8],
        (225000,
        2265.8610271903326,
        58.76602622336556,
        110834.05,
        0.21529937662783125,
        6.168374333872717,
        132315133.65660718,
        4.646199229999997) )
    np.testing.assert_almost_equal(
        tmp_res_Dic['tow_scenaria'][:8],
        (225000,
        2265.8610271903326,
        58.76602622336556,
        110834.05,
        0.21529937662783125,
        6.168374333872717,
        132315133.65660718,
        4.646199229999997) )
