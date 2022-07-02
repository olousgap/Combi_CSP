from CombiCSP import OutputContainer
import CombiCSP.SolarGeometry as sgh
import CombiCSP.CSP as cspC


class SolarTroughCalcs():
    def __init__(self,
        foc_len = 0.88 
        ,N = 1800 
        ,L = 25 
        ,Ws = 18 
        ,Wr = 0.07 
        ,Wc = 5.76 
        ):
        """_summary_

        Args:
            foc_len (float): [m] focal length CSPP T.1 in Mosleh19. Defaults to 0.88.
            N (int):  # [m * troughs] 25 * 48 CSPP pp.4 in Mosleh19 for 250 kWe turbine. Defaults to 1800.
            L (int): [m * troughs] 12 * 40 DISS pp.3 in Zarza04 for 70MWe turbine  . Defaults to 25.
            Ws (int): [m] width between rows 18 INDITEP in pp.6 Fraidenraich13,  pp.5 Zarza06. Defaults to 18.
            Wr (float): tube outer diameter [m]. Defaults to 0.07.
            Wc (float): collector width [m] 5.76 DISS pp.3 in Zarza04, 3.1 CSPP T.1 in Mosleh19 5-7.5 in SAM. Defaults to 5.76.
        """        
        

        self.foc_len = foc_len
        self.N = N
        self.L = L
        self.Ws = Ws
        self.Wr = Wr
        self.Wc = Wc 
        
    @property
    def area(self):
        return cspC.Ac(self.Wc, self.L, self.N)

    @property
    def Cg(self):
        return cspC.Cg_tro(Wc=self.Wc, Wr= self.Wr, L=self.L, N= self.N)

    def perform_calcs_EW(self, Ib, Tr=318, hoy=sgh.HOYS_DEFAULT):
        """Calculation for a solar trough oriented EW for a year per hour 

        Args:
            Ib (pd.Series): beam irradiance
            Tr (float, optional): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04. Defaults to 318.
            hoy (np.array, optional): _description_. Defaults to sgh.HOYS_DEFAULT.

        Returns:
            OutputContainer: Object that contains the power [MW] for each hour for the trough.
        """ 
        IAM=cspC.IAM_tro(hoy)
        data = cspC.di_sst(Ib=Ib,costhetai= cspC.costhetai_EW(),IAM=IAM, Tr=Tr, Wc=self.Wc, Wr=self.Wr, Ws=self.Ws, L=self.L, N=self.N)
        return OutputContainer(data = data, A_helio=self.area, Ctow=self.Cg)

    def perform_calcs_NS(self, Ib, Tr=318., hoy=sgh.HOYS_DEFAULT):
        """Calculation for a solar trough oriented NS for a year per hour 

        Args:
            Ib (pd.Series): beam irradiance
            Tr (float, optional): [oC] the working fluid temperature in the receiver, 350oC at DISS pp.3,7 in Zarza04. Defaults to 318.
            hoy (np.array, optional): _description_. Defaults to sgh.HOYS_DEFAULT.

        Returns:
            OutputContainer: Object that contains the power [MW] for each hour for the trough.
        """        
        IAM=cspC.IAM_tro(hoy)
        data = cspC.di_sst(Ib=Ib,costhetai=cspC.costhetai_NS(),IAM=IAM, 
                      Tr=Tr, 
                      Wc=self.Wc, Wr=self.Wr, Ws=self.Ws, 
                      L=self.L, N=self.N)
        return OutputContainer(data = data, A_helio=self.area, Ctow=self.Cg)
