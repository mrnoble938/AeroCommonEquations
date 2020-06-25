# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:41:29 2020

@author: noble

Aerospace Common Equations (ACE)
"""
import math

class Constants():
    def p0_psia() -> float:
        '''
        Sea Level Standard Day Pressure
        '''
        return 14.695949
    def t0_K() -> float:
        '''
        Sea Level Standard Day Temperature (K)
        '''
        return 288.15
    def rho0_slugperft3() -> float:
        '''
        Sea Level Standard Day Density
        '''
        return 0.0023769
    def a0_kts() -> float:
        '''
        Sea Level Standard Day Speed of Sound
        '''
        return 661.478827
    def g0_mps2() -> float:
        '''
        Sea Level Standard Gravitational Accel
        '''
        return 9.80665
    def meterperft() -> float:
        '''
        Meters per Foot
        '''
        return 0.3048
    def g0_fps2() -> float:
        '''
        Sea Level Standard Gravitational Accel
        '''
        return Constants.g0_mps2() / Constants.meterperft()
    def stdlapse_degCperft() -> float:
        '''
        Standard Lapse Rate
        '''
        return 0.0065 * Constants.meterperft()
    def stdlapse_degFperft() -> float:
        '''
        Standard Lapse Rate
        '''
        return Constants.stdlapse_degCperft() * (9/5)
    def tropopausealt_ft() -> float:
        '''
        Tropopause Altitude (ft)
        '''
        return 11000 / Constants.meterperft()
    def ftpermeter() -> float:
        '''
        Feet per meter
        '''
        return 1 / Constants.meterperft()
    def meterperinch() -> float:
        '''
        Meters per inch
        '''
        return Constants.meterperft() / 12
    def meterperNm() -> int:
        '''
        Meters per Nautical Mile
        '''
        return 1852
    def ftperNm() -> float:
        '''
        Feet per Nautical Mile
        '''
        return Constants.meterperNm() / Constants.meterperft()
    def Nmperft() -> float:
        '''
        Nautical Miles per foot
        '''
        return 1 / Constants.ftperNm()
    def ftperStMi() -> int:
        '''
        Feet per Statute Mile
        '''
        return 5280
    def StMiperft() -> float:
        '''
        Statute Miles per foot
        '''
        return 1 / Constants.ftperStMi()
    def kgperlb() -> float:
        '''
        Kilograms per pound
        '''
        return 0.45359237
    def lbperkg() -> float:
        '''
        Pounds per kilogram
        '''
        return 1 / Constants.kgperlb()
    def psiperPa() -> float:
        '''
        Pounds per Sq Inch per Pascale
        '''
        return Constants.meterperinch() ** 2 / \
            (Constants.g0_mps2() * Constants.kgperlb())
    def Paperpsi() -> float:
        '''
        Pascales per Pounds per Sq Inch
        '''
        return 1 / Constants.psiperPa()
    def PaperMillibar() -> int:
        '''
        Pascales per Millibar
        '''
        return 100
    def MillibarperPa() -> float:
        '''
        Millibar per Pascale
        '''
        return 1 / Constants.PaperMillibar()
    def PaperinHg() -> float:
        '''
        Pascales per inHg
        '''
        return 3386.389
    def inHgperPa() -> float:
        '''
        inHg per Pascale
        '''
        return 1 / Constants.PaperinHg()
    def fpmperkt() -> float:
        '''
        Feet per minute per knot
        '''
        return Constants.ftperNm() / 60
    def ktperfpm() -> float:
        '''
        Knots per feet per minute
        '''
        return 1 / Constants.fpmperkt()
    def fpsperkt() -> float:
        '''
        Feet per second per knot
        '''
        return Constants.fpmperkt() / 60
    def ktperfps() -> float:
        '''
        Knots per Feet per second
        '''
        return 1 / Constants.fpsperkt()
    def mpsperkt() -> float:
        '''
        Meters per second per knot
        '''
        return Constants.meterperNm() / 3600
    def ktpermps() -> float:
        '''
        Knots per Meters per Second
        '''
        return 1 / Constants.mpsperkt()
    def eps() -> float:
        '''
        A very small number
        '''
        return 0.0000001
    def sma_wgs84() -> float:
        '''
        Semi Major Axis of the earth in feet. Based on WGS84
        '''
        return 6378137 / Constants.meterperft()
    def flatfactor_wgs84() -> float:
        '''
        Flattening Factor describing ellipsoid of the earth. Based on WGS84
        '''
        return 1 / 298.257223563
    def b_wgs84() -> float:
        '''
        semi minor axis of the earth. Based on WGS84
        '''
        return Constants.sma_wgs84() * (1 - Constants.flatfactor_wgs84())
    def e2_wgs84() -> float:
        '''
        Eccentricity squared. Based on WGS84
        '''
        return Constants.flatfactor_wgs84() * \
            (2 - Constants.flatfactor_wgs84())
    
class AtmosphericProperties():
    def delta_atmo(hp_ft:float):
        '''
        Parameters
        ----------
        hp_ft : float
            Pressure altitude in Feet

        Returns
        -------
        Pressure Ratio
            delta.
        '''
        if hp_ft < Constants.tropopausealt_ft():
            delta = (1 - 6.87559e-06 * hp_ft) ** 5.25588
        else:
            delta = 0.22336 * \
                math.exp(-4.80634e-05 * (hp_ft - Constants.tropopausealt_ft()))
        return delta
    def pressure_alt(delta:float):
        '''
        Parameters
        ----------
        delta : float
            Pressure Ratio.

        Returns
        -------
        hp_ft: float
            pressure alt in feet
        '''
        if delta > \
            AtmosphericProperties.delta_atmo(Constants.tropopausealt_ft()):
            hp_ft = (1 - delta ** (1/5.25588)) / 6.87559e-06
        else:
            hp_ft = Constants.tropopausealt_ft() - (math.log(delta/0.22336) / \
                4.80634e-05)
        return hp_ft
    def delta_total(hp_ft:float,Mach:float):
        '''
        Parameters
        ----------
        hp_ft : float
            pressure alt in ft.
        Mach : float
            Mach number.

        Returns
        -------
        total pressure ratio.
        '''
        delta = AtmosphericProperties.delta_atmo(hp_ft)
        delta_total = delta * (1 + 0.2*Mach**2)**3.5
        return delta_total
    def hpfromPs_ft(staticPs:float):
        '''
        Parameters
        ----------
        staticPs : float
            Static Pressure in psia.

        Returns
        -------
        hp_ft : float
            pressure alt in ft.
        '''
        delta = staticPs / Constants.p0_psia()
        hp_ft = AtmosphericProperties.pressure_alt(delta)
        return hp_ft
    def Psfromhp_psia(hp_ft:float):
        '''
        Parameters
        ----------
        hp_ft : float
            Pressure Alt in ft.

        Returns
        -------
        staticPs : float
            Freestream static pressure.
        '''
        delta = AtmosphericProperties.delta_atmo(hp_ft)
        staticPs = delta * Constants.p0_psia()
        return staticPs
    def temp_ISA_C(hp_ft:float):
        '''
        Parameters
        ----------
        hp_ft : float
            Pressure Alt in ft.
        Returns
        -------
        temp_ISA : float
            ISA Temp in degC.
        '''
        if hp_ft < Constants.tropopausealt_ft():
            temp_ISA = 15 - Constants.stdlapse_degCperft() * hp_ft
        else:
            temp_ISA = -56.5
        return temp_ISA
    def temp_ISA_F(hp_ft:float):
        '''
        Parameters
        ----------
        hp_ft : float
            Pressure Alt in ft.
        Returns
        -------
        temp_ISA : float
            ISA temp in defF.
        '''
        if hp_ft < Constants.tropopausealt_ft():
            temp_ISA = 59 - Constants.stdlapse_degFperft() * hp_ft
        else:
            temp_ISA = -69.7
        return temp_ISA
    def dISA_C(OAT:float,hp_ft:float):
        '''
        Parameters
        ----------
        OAT : float
            Oatside air temp in degC.
        hp_ft : float
            Pressure Alt in ft.
        Returns
        -------
        dISA : float
            deviation from standard day temp degC.
        '''
        dISA = OAT - AtmosphericProperties.temp_ISA_C(hp_ft)
        return dISA
    def OAT_C(dISA:float,hp_ft:float):
        '''
        Parameters
        ----------
        dISA : float
            deviation from standard day temp degC.
        hp_ft : float
            Pressure Alt in ft.

        Returns
        -------
        OAT : float
            outside air temp in degC.
        '''
        if hp_ft < Constants.tropopausealt_ft():
            OAT = dISA + 15 - Constants.stdlapse_degCperft() * hp_ft
        else:
            OAT = dISA - 56.5
        return OAT
    def OAT_from_TAT_C(TAT:float,Mach:float,RF:float=1):
        '''
        Parameters
        ----------
        TAT : float
            total air temp in degC.
        Mach : float
            mach number.
        RF : TYPE, optional
            Recovery Factor, temp RF (less than 1) could be used
            to correct for flow not being perfectly adiabatic. 
            The default is 1.
        Returns
        -------
        OAT : TYPE
            outside air temp degC.
        '''
        RamRiseFactor = 1 + (0.2 * (Mach ** 2) * RF)
        OAT = ((TAT + 273.15) / RamRiseFactor) - 273.15
        return OAT
    def TAT_from_OAT_C(OAT:float,Mach:float,RF:float=1):
        '''
        Parameters
        ----------
        OAT : float
            outside air temp degC.
        Mach : float
            mach number.
        RF : float, optional
            Recovery Factor, temp RF (less than 1) could be used
            to correct for flow not being perfectly adiabatic. 
            The default is 1.
        Returns
        -------
        TAT : float
            Total Air Temp in degC.
        '''
        RamRiseFactor = 1 + (0.2 * (Mach ** 2) * RF)
        TAT = (OAT + 273.15) * RamRiseFactor - 273.15
        return TAT
    def theta_atmo(tempC:float):
        '''
        Parameters
        ----------
        tempC : float
            temperature in degC.
        Returns
        -------
        theta : float
            temperature ratio.
        '''
        theta = (273.15 + tempC) / Constants.t0_K()
        return theta
    def sigma_atmo(hp_ft:float,tempC:float):
        '''
        Parameters
        ----------
        hp_ft : float
            Pressure alt in ft.
        tempC : float
            temperature in degC.
        Returns
        -------
        sigma : float
            density ratio.
        '''
        sigma = AtmosphericProperties.delta_atmo(hp_ft) / \
            AtmosphericProperties.theta_atmo(tempC)
        return sigma
    def rho_slugperft3(hp_ft:float,OAT_C:float):
        '''
        Parameters
        ----------
        hp_ft : float
            pressure alt in ft.
        OAT_C : float
            oatside air temp in degC.
        Returns
        -------
        rho : float
            density of air [slug/ft^3].
        '''
        sigma = AtmosphericProperties.sigma_atmo(hp_ft, OAT_C)
        rho = Constants.rho0_slugperft3() * sigma
        return rho
    def accel_factor(Mach:float,hp_ft:float,ConstSpeedCode:int,OAT:float=None):
        '''
        Parameters
        ----------
        Mach : float
            mach number.
        hp_ft : float
            pressure alt in feet.
        ConstSpeedCode : int
            0 = Constant KCAS
            1 = Constant Mach
            2 = Constant KEAS
            3 = Constant KTAS.
        OAT : float, optional
            outside air temp in degC. The default is None.
        Returns
        -------
        AF : float
            Acceleration Factor.
        '''
        PHI = (((1 + 0.2 * (Mach ** 2)) ** 3.5) - 1) / \
            (0.7 * (Mach ** 2) * (1 + 0.2 * (Mach **2)) ** 2.5)
        if OAT == None:
            TempCorr = 1
        else:
            temp_ISA = AtmosphericProperties.temp_ISA_C(hp_ft)
            TempCorr = (273.15 + temp_ISA) / (273.15 + OAT)
        if ConstSpeedCode == 0:
            if hp_ft < Constants.tropopausealt_ft():
                AF = 1 + (0.7 * (Mach ** 2) * (PHI - 0.190263*TempCorr))
            else:
                AF = 1 + (0.7 * (Mach ** 2) * PHI)
        elif ConstSpeedCode == 1:
            if hp_ft < Constants.tropopausealt_ft():
                AF = 1 - (0.133184 * (Mach ** 2) * TempCorr)
            else:
                AF = 1
        elif ConstSpeedCode == 2:
            if hp_ft < Constants.tropopausealt_ft():
                AF = 1 + (0.7 * (Mach ** 2) * (1 - 0.190263 * TempCorr))
            else:
                AF = 1 + (0.7 * (Mach ** 2))
        elif ConstSpeedCode == 3:
            AF = 1
        else:
            print("ConstSpeedCode must be 0,1,2, or 3")
        return AF
    def dH_geo_ft(OAT:float,hp_ft:float,dHp_ft:float):
        '''
        Parameters
        ----------
        OAT : float
            oatside air temp [degC].
        hp_ft : float
            pressure altitude [ft].
        dHp_ft : float
            Pressure alt increment [ft].
        Returns
        -------
        dHGEO : float
            Geopotential Altitude Increment [ft].

        '''
        tempISA = AtmosphericProperties.temp_ISA_C(hp_ft)
        dHGEO = dHp_ft * ((273.15 + OAT) / (273.15 + tempISA))
        return dHGEO
    def dhp_ft(OAT:float,hp_ft:float,dHgeo_ft:float):
        '''
        Parameters
        ----------
        OAT : float
            oatside air temp [degC].
        hp_ft : float
            pressure alt [ft].
        dHgeo_ft : float
            geopotential alt increment [ft].
        Returns
        -------
        dHP : float
            pressure alt increment [ft].
        '''
        tempISA = AtmosphericProperties.temp_ISA_C(hp_ft)
        dHP = dHgeo_ft * ((273.15 + tempISA) / (273.15 + OAT))
        return dHP
    def dynamic_viscosity(OAT:float):
        '''
        Parameters
        ----------
        OAT : float
            outside air temp [degC].
        Returns
        -------
        mu : float
            dynamic viscosity of air [psi-sec].
        '''
        tempK = OAT + 273.15
        mu = (0.000001458 * (tempK ** 1.5)) / (tempK + 110.4) * \
            ((Constants.meterperinch() ** 2) / \
            (Constants.g0_mps2() * Constants.kgperlb()))
        return mu

class AirspeedConversions():
    def keas_from_kcas(kcas:float,hp_ft:float):
        '''
        Parameters
        ----------
        kcas : float
            calibrated airspeed [kt].
        hp_ft : float
            pressure altitude [ft].
        Returns
        -------
        keas : float
            equivalent airspeed [kt].
        '''
        delta = AtmosphericProperties.delta_atmo(hp_ft)
        keas = math.sqrt(5) * Constants.a0_kts() * \
            math.sqrt(((((((((1 + (0.2 * ((kcas / Constants.a0_kts()) ** 2))) \
                             ** 3.5) - 1) / \
                delta) + 1) ** (1/3.5)) - 1) * delta))
        return keas
    def kcas_from_keas(keas:float,hp_ft:float):
        '''
        Parameters
        ----------
        keas : float
            equivalent airspeed [kt].
        hp_ft : float
            pressure altitude [ft].
        Returns
        -------
        kcas : float
            calibrated airspeed [kt].
        '''
        delta = AtmosphericProperties.delta_atmo(hp_ft)
        kcas = math.sqrt((((((((((keas / (math.sqrt(5) * \
                    Constants.a0_kts())) ** 2) / delta) + 1) ** 3.5) \
                - 1) * delta + 1) ** (1/3.5)) - 1) / 0.2) * Constants.a0_kts()
        return kcas
    def ktas_from_keas(keas:float,hp_ft:float,OAT_C:float):
        '''        
        Parameters
        ----------
        keas : float
            equivalent airspeed [kt].
        hp_ft : float
            pressure altitude [ft].
        OAT_C : float
            outside air temp [degC].
        Returns
        -------
        ktas : float
            true airspeed [kt].
        '''
        sigma = AtmosphericProperties.sigma_atmo(hp_ft,OAT_C)
        ktas = keas / math.sqrt(sigma)
        return ktas
    def keas_from_ktas(ktas:float,hp_ft:float,OAT_C:float):
        '''
        Parameters
        ----------
        ktas : float
            true airspeed [kt].
        hp_ft : float
            pressure alt [ft].
        OAT_C : float
            outside air temp [degC].
        Returns
        -------
        keas : float
            equivalent airspeed [kt].
        '''
        sigma = AtmosphericProperties.sigma_atmo(hp_ft, OAT_C)
        keas = ktas * math.sqrt(sigma)
        return keas
    def ktas_from_kcas(kcas:float,hp_ft:float,OAT_C:float):
        '''
        Parameters
        ----------
        kcas : float
            calibrated airspeed [kt].
        hp_ft : float
            pressure altitude [ft].
        OAT_C : float
            outside air temp [degC].
        Returns
        -------
        ktas : float
            true airspeed [kt].
        '''
        keas = AirspeedConversions.keas_from_kcas(kcas, hp_ft)
        ktas = AirspeedConversions.ktas_from_keas(keas, hp_ft, OAT_C)
        return ktas
    def kcas_from_ktas(ktas:float,hp_ft:float,OAT_C:float):
        '''
        Parameters
        ----------
        ktas : float
            true airspeed [kt].
        hp_ft : float
            pressure alt [ft].
        OAT_C : float
            outside air temp [degC].
        Returns
        -------
        kcas : float
            calibrated airspeed [kt].
        '''
        keas = AirspeedConversions.keas_from_ktas(ktas, hp_ft, OAT_C)
        kcas = AirspeedConversions.kcas_from_keas(keas, hp_ft)
        return kcas
    def mach_from_deltap(pt:float,ps:float):
        '''
        Parameters
        ----------
        pt : float
            total pressure [psia].
        ps : float
            static pressure [psia].
        Returns
        -------
        mach : float
            mach number.
        '''
        mach = math.sqrt((5 * ((pt / ps) ** (1/3.5))) - 5)
        return mach
    def mach_from_keas(keas:float,hp_ft:float):
        '''
        Parameters
        ----------
        keas : float
            equivalent airspeed [kt].
        hp_ft : float
            pressure altitude [ft].
        Returns
        -------
        mach : float
            mach number.
        '''
        delta = AtmosphericProperties.delta_atmo(hp_ft)
        mach = keas / (Constants.a0_kts() * math.sqrt(delta))
        return mach
    def mach_from_kcas(kcas:float,hp_ft:float):
        '''
        Parameters
        ----------
        kcas : float
            calibrated airspeed [kt].
        hp_ft : float
            pressure alt [ft].
        Returns
        -------
        mach : float
            mach number.
        '''
        keas = AirspeedConversions.keas_from_kcas(kcas, hp_ft)
        mach = AirspeedConversions.mach_from_keas(keas, hp_ft)
        return mach
    def mach_from_ktas(ktas:float,OAT_C:float):
        '''
        Parameters
        ----------
        ktas : float
            true airspeed [kt].
        OAT_C : float
            outside air temp [degC].
        Returns
        -------
        mach : float
            mach number.
        '''
        mach = ktas / (Constants.a0_kts() * math.sqrt((OAT_C + 273.15) / \
                                                      288.15))
        return mach
    def keas_from_mach(mach:float,hp_ft:float):
        '''
        Parameters
        ----------
        mach : float
            mach number.
        hp_ft : float
            pressure alt [ft].
        Returns
        -------
        keas : float
            equivalent airspeed [kt].
        '''
        delta = AtmosphericProperties.delta_atmo(hp_ft)
        keas = mach * Constants.a0_kts() * math.sqrt(delta)
        return keas
    def kcas_from_mach(mach:float,hp_ft:float):
        '''
        Parameters
        ----------
        mach : float
            mach number.
        hp_ft : float
            pressure alt [ft].
        Returns
        -------
        kcas : float
            calibrated airspeed [kt].
        '''
        keas = AirspeedConversions.keas_from_mach(mach, hp_ft)
        kcas = AirspeedConversions.kcas_from_keas(keas, hp_ft)
        return kcas
    def keas_from_qpsf(q:float):
        '''
        Parameters
        ----------
        q : float
            dynamic pressure [psf].
        Returns
        -------
        keas : float
            equivalent airspeed [kt].
        '''
        keas = math.sqrt(295.375 * q)
        return keas
    def ktas_from_mach(mach:float,OAT_C:float):
        '''
        Parameters
        ----------
        mach : float
            mach number.
        OAT_C : float
            outside air temp [degC].
        Returns
        -------
        ktas : float
            true airspeed [kt].
        '''
        ktas = mach * Constants.a0_kts() * math.sqrt((OAT_C + 273.15) / 288.15)
        return ktas
    def qpsf_from_keas(keas:float):
        '''
        Parameters
        ----------
        keas : float
            equivalent airspeed [kt].
        Returns
        -------
        qpsf : float
            dynamic pressure [psf].
        '''
        qpsf = keas ** 2 / 295.375
        return qpsf
    def qpsf_from_mach(mach:float,hp_ft:float):
        '''
        Parameters
        ----------
        mach : float
            mach number.
        hp_ft : float
            pressure alt [ft].
        Returns
        -------
        qpsf : float
            dynamic pressure [psf].
        '''
        delta = AtmosphericProperties.delta_atmo(hp_ft)
        qpsf = 1481.352 * (mach ** 2) * delta
        return qpsf
    def kcas_from_qc(qc:float):
        '''
        Converts impact pressure in psid to kcas
        Parameters
        ----------
        qc : float
            impact pressure pt - ps (psid).
        Returns
        -------
        kcas : float
            calibrated airspeed [kt].
        '''
        kcas = math.sqrt(5) * Constants.a0_kts() * \
            math.sqrt((((qc / Constants.p0_psia()) + 1) ** (1/3.5)) - 1)
        return kcas
    def qc_from_kcas_psid(kcas:float):
        '''
        Parameters
        ----------
        kcas : float
            calibrated airspeed [kt].
        Returns
        -------
        qc : float
            impact pressure pt - ps (psid).
        '''
        qc = (((((kcas / math.sqrt(5) * Constants.a0_kts()) ** 2) + 1) ** 3.5) \
            - 1) * Constants.p0_psia()
        return qc
    def mps_from_kts(v_kts:float):
        '''
        Parameters
        ----------
        v_kt : float
            speed in knots.
        Returns
        -------
        v_mps : float
            speed in meters per sec.
        '''
        v_mps = Constants.mpsperkt() * v_kts
        return v_mps
    def kts_from_mps(v_mps:float):
        '''
        Parameters
        ----------
        v_mps : float
            speed in meters per sec.
        Returns
        -------
        v_kts : float
            speed in knots.
        '''
        v_kts = Constants.ktpermps() * v_mps
        return v_kts
    def fps_from_kts(v_kts:float):
        '''
        Parameters
        ----------
        v_kts : float
            speed in knots.
        Returns
        -------
        v_fps : float
            speed in feet per sec.
        '''
        v_fps = Constants.fpsperkt() * v_kts
        return v_fps
    def kts_from_fps(v_fps:float):
        '''
        Parameters
        ----------
        v_fps : float
            speed in feet per sec.
        Returns
        -------
        v_kts : float
            speed in knots.
        '''
        v_kts = Constants.ktperfps() * v_fps
        return v_kts
    def fpm_from_kts(v_kts:float):
        '''
        Parameters
        ----------
        v_kts : float
            speed in knots.
        Returns
        -------
        v_fpm : float
            speed in feet per minute.
        '''
        v_fpm = Constants.fpmperkt() * v_kts
        return v_fpm
    def kts_from_fpm(v_fpm:float):
        '''
        Parameters
        ----------
        v_fpm : float
            speed in feet per min.
        Returns
        -------
        v_kts : float
            speed in knots.
        '''
        v_kts = Constants.ktperfpm() * v_fpm
        return v_kts
    def mph_from_kts(v_kts:float):
        '''
        Parameters
        ----------
        v_kts : float
            speed in knots.
        Returns
        -------
        v_mph : float
            speed in mph.
        '''
        v_mph = (Constants.ftperNm() / Constants.ftperStMi()) * v_kts
        return v_mph
    def kts_from_mph(v_mph:float):
        '''
        Parameters
        ----------
        v_mph : float
            speed in miles per hour.
        Returns
        -------
        v_kts : float
            speed in knots.

        '''
        v_kts = (Constants.ftperStMi() / Constants.ftperNm()) * v_mph
        return v_kts

class UnitConversions():
    def meters_from_feet(dist_ft:float):
        dist_m = dist_ft * Constants.meterperft()
        return dist_m
    def feet_from_meters(dist_m:float):
        dist_ft = dist_m * Constants.ftpermeter()
        return dist_ft
    def feet_from_Nm(dist_Nm:float):
        dist_ft = dist_Nm * Constants.ftperNm()
        return dist_ft
    def feet_from_stmi(dist_stmi:float):
        dist_ft = dist_stmi * Constants.ftperStMi()
        return dist_ft
    def nm_from_stmi(dist_stmi:float):
        dist_nm = dist_stmi * Constants.ftperStMi() * Constants.Nmperft()
        return dist_nm