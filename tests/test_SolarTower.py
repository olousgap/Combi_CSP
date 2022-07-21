# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
"""
import pathlib
import pytest 
import numpy as np
import pandas as pd 

# import CombiCSP.SolarGeometry as sgh
# import CombiCSP.Transmittance as cspTr
# import CombiCSP.CSP as cspC
from CombiCSP import SolarTowerCalcs, SolarSystemLocation

@pytest.fixture
def st():
    """SolarTower ExampleData
    """    
    slobj = SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)
    return SolarTowerCalcs(alt = 200*10e-3 , Ht = 0.1, Ar = 99.3 , A_helio = 225000, slobj=slobj)

@pytest.fixture
def st2():
    """SolarTower ExampleData
    """    
    slobj = SolarSystemLocation(lat=12, lon=24, mer=-25, dt_gmt=+2, alt=0)
    return SolarTowerCalcs(alt = 200*10e-3 , Ht = 0.1, Ar = 99.3 , A_helio = 225000, slobj=slobj)


@pytest.fixture
def Ib():
    """SolarTower ExampleData
    """
    # pytest is configured in vscode at the root
    FNAME = pathlib.Path('tests/example_data/tmy_35.015_25.755_2005_2020.csv')
    pvgis = pd.read_csv(FNAME, header=16, nrows=8776-16, parse_dates=['time(UTC)'], engine='python') 
    Ib = pvgis.loc[:,'Gb(n)']
    return Ib


def test_perform_calcs_site1(st, Ib):
    oTow = st.perform_calc(Ib)
    vals = oTow.data.values
    assert vals.sum() == pytest.approx(110832.60989446181,1e-3) , 'Sum of values deos not match'
    assert np.abs(vals).sum() == pytest.approx(125038.6806487658 ,1e-3), 'Absolute sum failed'
    assert vals.std() == pytest.approx(18.822240024085243,1e-6) , 'Std is incorrect'
    assert np.abs(vals).std() == pytest.approx(17.6241589171352,1e-6) , 'Std is incorrect'
    
    # assert oTow.PowerMax_MW == pytest.approx(17.6241589171352,1e-6) , 'Max Power is incorrect'
    # assert oTow.PowerMax_CF == pytest.approx(17.6241589171352,1e-6) , 'CF is incorrect'
    # assert oTow.Energy_MWh == pytest.approx(17.6241589171352,1e-6) , 'Energy_MWh is incorrect'


def test_perform_calcs_st2(st2, Ib):
    """ This are calcualtaion in an alternate site (to see if the latitutes are used. )
    
    TODO This could be more precse

    Args:
        st2 (_type_): _description_
        Ib (_type_): _description_
    """    
    oTow = st2.perform_calc(Ib)
    vals = oTow.data.values
    assert not( vals.sum() == pytest.approx(110832.60989446181,1e-3)        ), 'Sum of values deos not match'
    assert not( np.abs(vals).sum() == pytest.approx(125038.6806487658 ,1e-3)) , 'Absolute sum failed'
    assert not( vals.std() == pytest.approx(18.822240024085243,1e-6)        ), 'Std is incorrect' 
    assert not( np.abs(vals).std() == pytest.approx(17.6241589171352,1e-6)  ), 'Std is incorrect' 
    
    # assert oTow.PowerMax_MW == pytest.approx(17.6241589171352,1e-6) , 'Max Power is incorrect'
    # assert oTow.PowerMax_CF == pytest.approx(17.6241589171352,1e-6) , 'CF is incorrect'
    # assert oTow.Energy_MWh == pytest.approx(17.6241589171352,1e-6) , 'Energy_MWh is incorrect'


