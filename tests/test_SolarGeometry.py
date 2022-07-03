# -*- coding: utf-8 -*- 
"""
    @Author: N. Papadakis
    @Date: 2022/07/02
"""

import pytest 
import numpy as np
import CombiCSP.SolarGeometry as sgh
from CombiCSP.solar_system_location import SolarSystemLocation

@pytest.fixture
def hoy_ex1():
    """hoy example data
    """    
    return np.array([1])
@pytest.fixture
def Crete():
    """Example site
    """    
    return SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)

@pytest.fixture
def R():
    """R for transmittance functions 
    units are unknown. For verifications see : 
    H.C. Hottel, A simple model for estimating the transmittance of direct solar radiation through clear atmospheres, 
    Solar Energy. 18 (1976) 129–134."
    """
    return 1

    
     
class Test_AM:
    """testing the air mass functions
    """    
    def test_AM(self, hoy_ex1):
        
        assert sgh.AM(hoy_ex1)[0] == pytest.approx(  -1.0308066793553488, abs=1e-3)
        assert sgh.air_mass(hoy_ex1, method='wiki')[0] == pytest.approx(  -1.0308066793553488, abs=1e-3)
        assert sgh.air_mass(hoy_ex1 )[0] == pytest.approx(  -1.0308066793553488, abs=1e-3)
        # pytest.approx(sgh.AM(hoy_ex1, method='pveducation')[0], -13.37218176, abs=1e-3)
        
    def test_AM2(self,hoy_ex1):
        assert sgh.AM2(hoy_ex1)[0] == pytest.approx( -1.1149391977788454, abs=1e-3)
        assert sgh.air_mass(hoy_ex1, method='Kasten')[0] == pytest.approx( -1.1149391977788454, abs=1e-3)
        
    def test_AM3(self,hoy_ex1 ):
        assert sgh.AM3(hoy_ex1)[0] == pytest.approx( -1.1184591949819083, abs=1e-3)
        assert sgh.air_mass(hoy_ex1, method='Kasten-Young')[0] == pytest.approx( -1.1184591949819083, abs=1e-3)
        
    def test_AM4(self,hoy_ex1 ):
        assert sgh.AM4(hoy_ex1)[0] == pytest.approx( 1374.4966167568873, abs=1e-3)
        assert sgh.air_mass(hoy_ex1, method='Schoenberg')[0] == pytest.approx( 1374.4966167568873, abs=1e-3)
        


class Test_EOT_zen_ele:
    """testing the air mass functions
    """    
    def test_EOT(self, hoy_ex1, Crete):
        expected = -2.9041689600000002 
        assert sgh.EoT(hoy_ex1)[0] == pytest.approx(  expected , abs=1e-3)
         
    def test_tsol(self,hoy_ex1, Crete):
        expected = 0.551597184 
        assert sgh.tsol(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        assert Crete.tsol(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        
    def test_ele(self, hoy_ex1, Crete ):
        expected = -1.3257002183993918 
        assert sgh.ele(hoy_ex1)[0] == pytest.approx(expected, abs=1e-3)
        assert Crete.ele(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        
    def test_z(self,hoy_ex1, Crete ):
        expected =  2.8964965451942883 
        assert sgh.z(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        assert Crete.z(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)

    def test_azim(self,hoy_ex1, Crete):
        expected =  -0.577725020643127
        assert sgh.azim(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        assert Crete.azim(hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)

