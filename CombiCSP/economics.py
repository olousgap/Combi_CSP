# -*- coding: utf-8 -*- 
"""
    @Author: G.Arnaoutakis, N.Papadakis
    @Date: 2022/06/30
    @Credit: the functions were originally in  CSPEcon.py 
"""
#%%
import numpy as np
import pandas as pd
import numpy_financial as npf

from CombiCSP import OutputContainer

def discounted_payback_period(rate, cash_flows): # #https://sushanthukeri.wordpress.com/2017/03/29/discounted-payback-periods/
    '''Return the Discounted Payback Period'''
    cf_df = pd.DataFrame(cash_flows, columns=['UndiscountedCashFlows'])
    cf_df.index.name = 'Year'
    cf_df['DiscountedCashFlows'] = npf.pv(rate, pmt=0, nper=cf_df.index, fv=-cf_df['UndiscountedCashFlows'])
    cf_df['CumulativeDiscountedCashFlows'] = np.cumsum(cf_df['DiscountedCashFlows'])
    final_full_year = cf_df[cf_df.CumulativeDiscountedCashFlows < 0].index.values.max()
    fractional_yr = -cf_df.CumulativeDiscountedCashFlows[final_full_year ]/cf_df.DiscountedCashFlows[final_full_year + 1]
    payback_period = final_full_year + fractional_yr
    #bar(lifetime, cf_df['CumulativeDiscountedCashFlows'], color='c')
    #xlabel('Time (years)'), ylabel('€')
    #grid(linestyle="--", linewidth=0.5, color='.25', zorder=-10)
    #ylim([-1e8,1.25e8])
    return payback_period

def irr(values): # http://nullege.com/codes/search/numpy.irr
    """
    Return the Internal Rate of Return (IRR).
     
    This is the "average" periodically compounded rate of return
    that gives a net present value of 0.0; for a more complete explanation,
    see Notes below.
     
    Parameters
    ----------
    values : array_like, shape(N,)
        Input cash flows per time period.  By convention, net "deposits"
        are negative and net "withdrawals" are positive.  Thus, for example,
        at least the first element of `values`, which represents the initial
        investment, will typically be negative.
     
    Returns
    -------
    out : float
        Internal Rate of Return for periodic input values.
     
    Notes
    -----
    The IRR is perhaps best understood through an example (illustrated
    using np.irr in the Examples section below).  Suppose one invests
    100 units and then makes the following withdrawals at regular
    (fixed) intervals: 39, 59, 55, 20.  Assuming the ending value is 0,
    one's 100 unit investment yields 173 units; however, due to the
    combination of compounding and the periodic withdrawals, the
    "average" rate of return is neither simply 0.73/4 nor (1.73)^0.25-1.
    Rather, it is the solution (for :math:`r`) of the equation:
     
    .. math:: -100 + \\frac{39}{1+r} + \\frac{59}{(1+r)^2}
     + \\frac{55}{(1+r)^3} + \\frac{20}{(1+r)^4} = 0
     
    In general, for `values` :math:`= [v_0, v_1, ... v_M]`,
    irr is the solution of the equation: [G]_
     
    .. math:: \\sum_{t=0}^M{\\frac{v_t}{(1+irr)^{t}}} = 0
     
    References
    ----------
    .. [G] L. J. Gitman, "Principles of Managerial Finance, Brief," 3rd ed.,
       Addison-Wesley, 2003, pg. 348.
     
    Examples
    --------
    >>> round(irr([-100, 39, 59, 55, 20]), 5)
    0.28095
    >>> round(irr([-100, 0, 0, 74]), 5)
    -0.0955
    >>> round(irr([-100, 100, 0, -7]), 5)
    -0.0833
    >>> round(irr([-100, 100, 0, 7]), 5)
    0.06206
    >>> round(irr([-5, 10.5, 1, -8, 1]), 5)
    0.0886
     
    (Compare with the Example given for numpy.lib.financial.npv)
     
    """
    res = np.roots(values[::-1])
    mask = (res.imag == 0) & (res.real > 0)
    if res.size == 0:
        return np.nan
    res = res[mask].real
    # NPV(rate) = 0 can have more than one solution so we return
    # only the solution closest to zero.
    rate = 1.0/res - 1
    rate = rate.item(np.argmin(np.abs(rate)))
    return rate

def cashflow(Ecsp,csp_price,fuel_energy,eff,fuel_price,capital): 
    return round(Ecsp*csp_price+fuel_energy/eff*fuel_price-0.04*capital,2)

def expenses(Eaux, eff,fuel_price,capital):
    """_summary_

    Args:
        Eaux (_type_): auxilliary energy required to operate the facility 
        eff (_type_): fuel efficiency
        fuel_price (_type_): fuel price  [?]
        capital (_type_): _description_

    Returns:
        _type_: expenses
    """     
    return round(Eaux/eff*fuel_price-0.04*capital,2)

# %%
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



class Economic_environment():
    """this is a class that is responsible for containing basic economic parameters for fiscal analysis
    
    
    """    
    def __init__(self, 
        oil_price:float, Eoil:float,
        currency_units:str ='USD'):
        """_summary_

        Args:
            oil_price (float): _description_
            Eoil (float): _description_
            currency_units (str, optional): _description_. Defaults to 'USD'.
        """        
        
        self._oil_price = oil_price
        self._Eoil = Eoil
        self.currency_units = currency_units

    @property
    def oil_price(self):
        """return oil price

        Returns:
            float: _description_
        """        
        return self._oil_price

    def MWh_to_BOE(self, energy_MWh):
        """converts MWh to Barrel Oil Equivalent

        # [BOE] 1MWh = 0.5883BOE https://www.convert-me.com/en/convert/energy/kwh/kwh-to-boe.html?u=kwh&v=1%2C000
        
        Args:
            energy_MWh (_type_): _description_

        Returns:
            _type_: _description_
        """
        return energy_MWh*0.588

    def MWh_to_mBTU(self, energy_MWh):
        """converts MWh to   Equivalent natural gas m^3

        [m BTU] 1MWh = 3.412mBTU https://www.convert-me.com/en/convert/energy/kwh/kwh-to-mymmbtu.html?u=kwh&v=1%2C000

        Args:
            energy_MWh (_type_): _description_

        Returns:
            float: 
        """
        return energy_MWh*3.412

    
    def economics_for_Solar_tower(self, 
            oTow:OutputContainer,
            csp_area_costs,
            csp_energy_price,
            csp_discount_rate,
            power_block_cost,
            capital_csp,
        lifetime=range(30)):
        """This function performs an economic analysis on the performance output of a csp

        Args:
            oTow (OutputContainer): _description_
            csp_area_costs (_type_): _description_
            csp_energy_price (_type_): _description_
            csp_discount_rate (_type_): _description_
            capital_csp (float): description
            power_block_cost (_type_): _description_
            Eoil (_type_): _description_
            oil_price (_type_): _description_
            lifetime (_type_, optional): _description_. Defaults to range(30).

        Returns:
            _type_: _description_
        """    

        capital_csp_tow = oTow.A_helio* csp_area_costs + oTow.PowerMax_MW*power_block_cost
        revenue_csp_tow = cashflow(oTow.Energy_MWh, 
            csp_energy_price, self._Eoil, 
            0.4,
            -self.oil_price, capital_csp_tow)
        cash_flow_tow = [-capital_csp_tow] + [revenue_csp_tow for i in lifetime]
        dpb_tow = discounted_payback_period(csp_discount_rate, cash_flow_tow)
        npv_csp_tow = npf.npv (csp_discount_rate, [-capital_csp_tow] + [cashflow(oTow.Energy_MWh,csp_energy_price,self._Eoil,0.4,-self.oil_price,capital_csp_tow) for i in lifetime])
        irr_csp_tow = npf.irr([-capital_csp] + [cashflow(oTow.Energy_MWh,csp_energy_price,self._Eoil,0.4,-self.oil_price,capital_csp_tow) for i in lifetime])
        return {
            'A_helio': oTow.A_helio,
            'cash_flow_tow':cash_flow_tow,
            'tow_scenaria': (oTow.A_helio, oTow.Ctow, 
                            oTow.PowerMax_MW, oTow.Energy_MWh,
                            oTow.CF, dpb_tow, 
                            npv_csp_tow, irr_csp_tow, 
                            cash_flow_tow)
        }
        
    def economics_for_SolarTrough(self, 
        oTr:OutputContainer,
        csp_area_costs,
        csp_energy_price,
        csp_discount_rate,
        power_block_cost,
        lifetime=range(30)):
        """This function performs an economic analysis on the performance output of a csp

        Args:
            oTow (OutputContainer): _description_
            csp_area_costs (_type_): _description_
            csp_energy_price (_type_): _description_
            csp_discount_rate (_type_): _description_
            capital_csp (float): description
            power_block_cost (_type_): _description_
            Eoil (_type_): _description_
            oil_price (_type_): _description_
            lifetime (_type_, optional): _description_. Defaults to range(30).

        Returns:
            _type_: _description_
        """    
        capital_csp_tro = oTr.A_helio* csp_area_costs + oTr.PowerMax_MW*power_block_cost
        revenue_csp_tro = cashflow(oTr.Energy_MWh,csp_energy_price,self._Eoil,0.4,-self.oil_price,capital_csp_tro)

        cash_flow_tro = [-capital_csp_tro] + [revenue_csp_tro for i in range(30)]
        dpb_tro = discounted_payback_period(csp_discount_rate, cash_flow_tro)
        npv_csp_tro = npf.npv(csp_discount_rate, [-capital_csp_tro] \
            + [revenue_csp_tro for i in range(30)])
        irr_csp_tro = npf.irr([-capital_csp_tro] \
            + [revenue_csp_tro for i in range(30)])
        return {
            'A_helio': oTr.A_helio,
            'cash_flow_tow':cash_flow_tro,
            'tow_scenaria': (oTr.A_helio, oTr.Ctow, 
                        oTr.PowerMax_MW, oTr.Energy_MWh,
                        oTr.CF, dpb_tro, 
                        npv_csp_tro, irr_csp_tro, 
                        cash_flow_tro)}
# %%
