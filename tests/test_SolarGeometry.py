import pytest 
import numpy as np
import CombiCSP.SolarGeometry as sgh

@pytest.fixture
def hoy_ex1():
    """hoy example data
    """    
    return np.array([1])
    
class Test_da:
    def test_d1(self):
        h = np.array([1])
        pytest.approx(sgh.d1(h)[0], -13.37218176, abs=1e-3)
        pytest.approx(sgh.eda(h, method='pveducation')[0], -13.37218176, abs=1e-3)
        
    def test_d2(self):
        h = np.array([1])
        pytest.approx(sgh.d2(h)[0], -13.37218176, abs=1e-3)
        pytest.approx(sgh.eda(h, method='Katsaprakakis')[0], -13.37218176, abs=1e-3)
        
    def test_d3(self):
        h = np.array([1])
        pytest.approx(sgh.d3(h)[0], -13.37218176, abs=1e-3)
        pytest.approx(sgh.eda(h, method='-81')[0], -13.37218176, abs=1e-3)
        
    def test_d(self):
        h = np.array([1])
        pytest.approx(sgh.d(h)[0], -13.37218176, abs=1e-3)
        pytest.approx(sgh.eda(h, method='wiki')[0], -13.37218176, abs=1e-3)
        pytest.approx(sgh.eda(h)[0], -13.37218176, abs=1e-3)
    
def test_eda_raise_ValueError(hoy_ex1):
    with pytest.raises(ValueError):
        sgh.eda(hoy_ex1,method='noncense')
    sgh.eda(hoy_ex1,method='wiki')
    sgh.eda(hoy_ex1,method='Katsaprakakis')
    sgh.eda(hoy_ex1,method='-81')
    sgh.eda(hoy_ex1,method='pveducation')
    