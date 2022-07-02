# -*- coding: utf-8 -*- 
"""
    @Author: G. Arnaoutakis
    @Date: 2022/06/30
"""
#%%
import numpy as np
import pandas as pd
import numpy_financial as npf

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
def expenses(eff,fuel_price,capital): return round(Eaux/eff*fuel_price-0.04*capital,2)

# %%
