import pytest 
import numpy as np
import CombiCSP.SolarGeometry as sgh

@pytest.fixture
def hoy_ex1():
    """hoy example data
    """    
    return np.array([1])
    
class Test_da:
    def test_d1(self, hoy_ex1):
        assert sgh.d1(hoy_ex1)[0] == pytest.approx(-13.37218176, abs=1e-3)
        assert sgh.eda(hoy_ex1, method='pveducation')[0] == pytest.approx(-13.37218176, abs=1e-3)
        # pytest.approx(sgh.d1(h)[0], -13.37218176, abs=1e-3)
        # pytest.approx(sgh.eda(h, method='pveducation')[0], -13.37218176, abs=1e-3)
        
    def test_d2(self, hoy_ex1):
        assert sgh.d2(hoy_ex1)[0] == pytest.approx( 23.41331062609803, abs=1e-3)
        assert sgh.eda(hoy_ex1, method='Katsaprakakis')[0] == pytest.approx( 23.41331062609803, abs=1e-3)
        
    def test_d3(self, hoy_ex1):
        assert sgh.d3(hoy_ex1)[0] == pytest.approx( 23.413310626097925, abs=1e-3)
        assert sgh.eda(hoy_ex1, method='-81')[0] == pytest.approx( 23.413310626097925, abs=1e-3)
        
    def test_d(self,hoy_ex1):
        assert sgh.d(hoy_ex1)[0] == pytest.approx( -0.4005513174270792, abs=1e-3)
        assert sgh.eda(hoy_ex1, method='wiki')[0] == pytest.approx( -0.4005513174270792, abs=1e-3)
        assert sgh.eda(hoy_ex1)[0] == pytest.approx( -0.4005513174270792, abs=1e-3)
    
def test_eda_raise_ValueError(hoy_ex1):
    with pytest.raises(ValueError):
        sgh.eda(hoy_ex1,method='noncense')
    sgh.eda(hoy_ex1,method='wiki')
    sgh.eda(hoy_ex1,method='Katsaprakakis')
    sgh.eda(hoy_ex1,method='-81')
    sgh.eda(hoy_ex1,method='pveducation')
    
    
     
class Test_AM:
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
        
        