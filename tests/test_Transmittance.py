import pytest 
import numpy as np
import CombiCSP.SolarGeometry as sgh
import CombiCSP.Transmittance as cspTr

@pytest.fixture
def hoy_ex1():
    """hoy example data
    """    
    return np.array([1])

@pytest.fixture
def Crete():
    """Example site
    """    
    return sgh.SolarSystemLocation(lat=35, lon=24, mer=-25, dt_gmt=+2, alt=0)

@pytest.fixture
def R():
    """R for transmittance functions 
    units are unknown. For verifications see : 
    H.C. Hottel, A simple model for estimating the transmittance of direct solar radiation through clear atmospheres, 
    Solar Energy. 18 (1976) 129–134."

    NP: I don't know that R is, how it is defined and what are its units
        this is used only as a unit test to make sure that the results dont change the `trig` functions to `np.trig` functions
    """
    return 1

@pytest.fixture
def alt():
    """
    NP: I suspect alt is the altitude of the system, but I don't know how to verify and what are its units
        this is used only as a unit test to make sure that the results dont change the `trig` functions to `np.trig` functions
    """    
    return 1
@pytest.fixture
def Ht():
    """ NP: I don't know that Ht is and what are its units
        this is used only as a unit test to make sure that the results dont change the `trig` functions to `np.trig` functions
    """    
    return 1


class Test_TransmittivityFunctions:
    """testing the Transmittivity functions

    **IMPORTANT**: These tests were developed only to make sure that the original implementation by George did not get
    affected by the change the `trig` functions to `np.trig` functions

    #TODO Everyone of this functions should be carefully validated at a later point.
    #TODO the Tr23m should have a different formula
    """    
    def test_Tr23km(self,  hoy_ex1, alt):

        expected = 0.4207775856677556
        assert cspTr.Tr23km(alt, hoy_ex1)[0] == pytest.approx(  expected , abs=1e-3)
         
    def test_Tr5km(self,hoy_ex1, alt):
        expected = 0.2930476046606646 
        assert cspTr.Tr5km(alt=alt, hoy=hoy_ex1)[0] == pytest.approx( expected, abs=1e-3)
        
        
    def test_TrD23km(self, R ):
        expected = 9.7184
        assert cspTr.TrD23km(R) == pytest.approx(expected, abs=1e-3)
        
        
    def test_TrD5km(self, R ):
        expected = 25.378999999999998
        assert cspTr.TrD5km(R) == pytest.approx( expected, abs=1e-3)
        

    def test_TrV23km(self, R):
        expected =  0.902815
        assert cspTr.TrV23km(R) == pytest.approx( expected, abs=1e-3)
        

    def test_TrV5km(self, R ):
        expected =  0.7462099999999999
        assert cspTr.TrV5km(R) == pytest.approx( expected, abs=1e-3)
        

    def test_TrVH(self,Ht, R, alt):
        expected =  0.9073711227199887
        assert cspTr.TrVH(Ht, R, alt) == pytest.approx( expected, abs=1e-3)
        