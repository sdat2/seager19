"""src.models.atmos.

This model was initially written in python by Dr. Naomi Henderson.

It includes the atmospheric Matsuno-Gill model, and the Seager1991 surface flux model.

It was refactored by Simon Thomas into a class structure / to be pylint compatible,
and to take the cfg DictConfig struct as input.

It includes both the tropical Atmospheric model
(Matsumo-Gill) with a single baroclinic mode,
and the surface fluxes calculated based on the Ocean model.

An introduction to the Matsumo-Gill model:

https://www.atmos.washington.edu/~dargan/591/591_8.pdf

pytest src/test/test_atmos.py
python3 src/models/atmos.py

Model Solution Method from Seager et al. 2019 Appendix:

The atmosphere equations are solved by Fourier transforming in longitude,
forming an equation for v for each zonal wavenumber that is finite
differenced, and the resulting tri-diagonal system is solved by matrix
inversion, transforming back into longitude. Finally, u and Φ are derived
by back-substitution. The ocean equations are solved using the ‘INC’
scheme, integrating the model forward, after spin-up with
climatological conditions, forced by the time-varying ECMWF wind stress
and, for the case with CO2 forcing, changing 𝑓′1 in the net surface
longwave radiation calculation. Change over 1958–2017 is computed by a
linear trend. The atmosphere model is solved forced by a Ts comprised of
the climatological mean for 1958–2017 plus and minus half of the SST trend
and the difference of the two simulations taken to derive the change.
For the coupled model, the ocean model is first forced with the change in
CO2 and climatological wind stress over 1958–2017. The resulting SST
trend, plus the imposed heating change over land, are used to force the
atmosphere model. The ocean model is forced again with both the changed
wind stress and the CO2 increase to derive a new SST change over 1958–2017
that is then used to force the atmosphere model. This iterative coupling
is repeated until equilibrium is reached, which takes just a few times.
There is a unique solution for any given value of CO2. The model wind
stress change is computed as 𝜌a𝑐D𝑊𝐮, where cD is a drag coefficient
and 𝐮 is the vector surface wind change computed by the atmosphere model,
which is added to the ECMWF climatological stresses. Since the atmosphere
model dynamics are only applicable in the tropics, the computed wind
stress anomaly is only applied to the ocean model between 20 S and 20 N,
and is linearly tapered to zero at 25 S and 2 N.

Example:
    Old code to record the initial parameter settings::

        # ------------- constants -----------------------
        # begining TCAM
        self.pr_max: float = 20.0 / 3600 / 24
        # 20 / seconds in hour / hours in day.
        self.relative_humidity: float = 0.80  # relative humidity uniformly 0.8
        self.number_iterations: int = 50  #  int
        self.R_earth.to_value() = 6.37e6  # metres
        self.sec_in_day = 86400  # seconds in day.
        self.omega_2 = 2 * (2 * np.pi / self.sec_in_day)  # 2 * rad per second
        self.latent_heat_vap = 2.5e6  # latent heat # J kg-1
        self.cp_air = 1000  #  cp_air is the specific heat capacity of air.
        # J kg-1 K-1
        self.b_coeff = (
            self.gravity * np.pi / (self.nbsq * self.theta_00 * self.height_tropopause)
        )
        self.eps_days = 0.75
        # Over 1958–2017, the CO2 changed from ~300 to ~400 ppm,
        # which would be about 0.75 W m−2
        # eps days might be the efficiency of entrainment.
        self.eps = 1.0 / (self.eps_days * self.sec_in_day)  # 1/.75 d
        self.eps_u = self.eps  # 1/.75 d
        self.eps_v = self.e_frac * self.eps  # e_frac=1/2 in paper
        # Newtonian cooling, K
        self.newtonian_cooling_coeff_k1 = self.b_coeff / (self.k_days * self.sec_in_day)
        self.eps_p = (np.pi / self.height_tropopause) ** 2 / (
            self.nbsq * self.k_days * self.sec_in_day
        )
        self.beta = self.omega_2 / self.R_earth.to_value()
        self.rho_air = 1.225  # kg m-3 - also called rho_00
        self.c_e = 0.00125  # 1.25e-3 # cE is an exchange coefficient
        self.emmisivity = 0.97  # problem here with second definintion.
        self.stefan_boltzman_const = 5.67e-8
        self.p_s = 1000  # pressure at the surface?
        self.es_0 = 6.11
        self.delta_temp = 1.0  # ΔT = 1 K
        self.f2 = 0.05  # f2 = 0.05
        # 'a_cloud_const' should decrease when deep convection happens above 28 degC
        # python main.py --multirun oc.nummode=6,7,8
        # a_cloud_const = Ts-temp_0_c;
        # a_cloud_const[a_cloud_const>28] = 40;
        # a_cloud_const[a_cloud_const<=28] = 80;
        # a_cloud_const = 0.01*a_cloud_const
        self.a_cloud_const = 0.6  # this isn't the option used in the paper.
        # basic parameters
        self.temp_0_c = 273.15  # zero degrees in kelvin
        self.f1_bar = 0.39  # f1 = 0.39
        # f'1  is the anomaly in f1—a parameter that can be adjusted
        # to control the variation in surface longwave radiation due
        # to a_cloud_const change in CO2
        self.u_bar = 5.0  # average velocity?
        self.temp_surface_bar: float = self.temp_0_c + 25  # 25C in Kelvin
        self.c_bar = 0.6  # C is the cloud cover. perhaps C_bar is the average.

"""
from typing import Tuple, Union, Any
import os
import numpy as np
from scipy.interpolate import interp2d
from scipy.fftpack import fft, ifft
from scipy.constants import zero_Celsius
from astropy.constants import R_earth
import matplotlib.pyplot as plt
import xarray as xr
from typeguard import typechecked
from omegaconf import DictConfig
from src.models.model_setup import ModelSetup
from src.utils import timeit
from src.constants import MODEL_NAMES, VAR_DICT


class Atmos:
    """Atmos class."""

    def __init__(self, cfg: DictConfig, setup: ModelSetup) -> None:
        """Initialise the atmos function

        Args:
            cfg (DictConfig): Config params to feed in.
            setup (ModelSetup): file structure for run.
        """
        self.atm = cfg.atm
        self.setup = setup
        self.it = 0

        # make axes
        self.x_axis = np.linspace(0, 360 - self.atm.dx, self.atm.nx)  # degrees
        self.y_axis_v = np.linspace(
            self.atm.y_south_lim + self.atm.dy / 2,
            self.atm.y_north_lim - self.atm.dy / 2,
            self.atm.ny,
        )  # degrees
        self.y_axis_u = np.linspace(
            self.atm.y_south_lim + self.atm.dy,
            self.atm.y_north_lim - self.atm.dy,
            self.atm.ny - 1,
        )  # degrees
        self.y_axis_i = np.linspace(
            self.atm.y_south_lim + 3 * self.atm.dy / 2,
            self.atm.y_north_lim - 3 * self.atm.dy / 2,
            self.atm.ny - 2,
        )  # degrees

        # adding coriolis params at different y axis locations
        self.fcu = self.f_cor(self.y_axis_u)  # vector of coriolis force coefficient.

        # properties derived from grid axes
        self.x_spacing = self.x_axis[1] - self.x_axis[0]  # degrees
        self.y_spacing = self.y_axis_v[1] - self.y_axis_v[0]  # degrees
        self.dxm = self.x_spacing * R_earth.to_value() * np.pi / 180
        self.dym = self.y_spacing * R_earth.to_value() * np.pi / 180
        self.dym_2 = self.dym * self.dym

        # need to have the correct ordering of the wave numbers for fft
        if self.atm.nx % 2 == 0:
            self.kk_wavenumber = np.asarray(
                list(range(0, self.atm.nx // 2))
                + [0]
                + list(range(-self.atm.nx // 2 + 1, 0)),
                np.float64,
            )
        else:
            self.kk_wavenumber = np.asarray(
                list(range(0, (self.atm.nx - 1) // 2))
                + [0]
                + list(range(-(self.atm.nx - 1) // 2, 0)),
                np.float64,
            )

        # the different model names in a dict? - used by key from self.mem.
        self.names: dict = MODEL_NAMES

        # dict of variables that are read in.
        self.var: dict = VAR_DICT
        # temperature of the surface, cloud area fraction, surface wind, rel humidity.

        # END INIT.

    @typechecked
    def f_cor(self, y_axis: np.ndarray) -> np.ndarray:
        """Corriolis force coeff. Makes beta plane approximation.

        omega_2 = 2 * (2 * np.pi / sec_in_day) # 2 * rad per second

        Args:
            y (np.ndarray): latitude

        Returns:
            np.ndarray: Corriolis force coeff.
        """
        return self.atm.omega_2 * y_axis * np.pi / 180  ## convert from degrees to rad.

    # --------------- fluxes -----------------------------

    @typechecked
    def f_es(self, temperature: xr.DataArray) -> xr.DataArray:
        """Flux es.

        Args:
            temperature (xr.DataArray): temp.

        Returns:
            xr.DataArray: Flux es.
        """
        return self.atm.es_0 * np.exp(
            17.67 * (temperature - zero_Celsius) / (temperature - zero_Celsius + 243.5)
        )

    @typechecked
    def f_qs(self, temperature: xr.DataArray) -> xr.DataArray:
        """Flux q_s.

        q_s(Ts) is the saturation-specific humidity at the SST
        rqs(Ts) = q_s(Ts)

        Args:
            temperature (xr.DataArray): temp.

        Returns:
            xr.DataArray: Flux q_s.
        """
        return self.atm.e_factor * self.f_es(temperature) / self.atm.p_s

    @typechecked
    def f_dqs_dtemp(self, temperature: xr.DataArray) -> xr.DataArray:
        """Flux dqs dtemp.

        change in saturation humdity of the suface by temperature.

        Args:
            temperature (xr.DataArray): temp.

        Returns:
            xr.DataArray: flux q_s.
        """
        return (
            self.f_qs(temperature)
            * (17.67 * 243.5)
            / (temperature - zero_Celsius + 243.5) ** 2
        )

    @typechecked
    def f_qlh(
        self,
        temperature: xr.DataArray,
        u_sp: Union[xr.DataArray, float],
        rh_loc: xr.DataArray,
    ) -> xr.DataArray:
        """Heat flux from latent heat.

        It is assumed that the surface heat flux anomaly is
        dominated by longwave and latent heat fluxes and
        that the solar radiation does not change and the
        sensible heat flux anomaly is small.

        Args:
            temperature (xr.DataArray): the temperature dataarray.
            u_sp (Union[xr.DataArray, float]): the wind speed.
            rh_loc (xr.DataArray): relative humidity

        Returns:
            xr.DataArray: flux qlh.

        """
        return self.atm.qlh_coeff * u_sp * self.f_qs(temperature) * (1 - rh_loc)

    @typechecked
    def f_dqlh_dtemp(
        self,
        temperature: xr.DataArray,
        u_sp: Union[xr.DataArray, float],
        rh_loc: xr.DataArray,
    ) -> xr.DataArray:
        """Flux dqlh_dtemp.

        Args:
            temperature (xr.DataArray): temperature (in kelvin).
            u_sp (Union[xr.DataArray, float]): u speed.
            rh_loc (xr.DataArray): relative humidity.

        Returns:
            xr.DataArray: flux dqlh_dtemp.
        """
        return self.atm.qlh_coeff * u_sp * self.f_dqs_dtemp(temperature) * (1 - rh_loc)

    @typechecked
    def f_temp_a(self, temperature: xr.DataArray) -> xr.DataArray:
        """Temperature anomaly.

        Delta temp = 1 in paper.

        Args:
            temperature (xr.DataArray): temperature in kelvin.

        Returns:
            xr.DataArray: temperature anomaly.
        """
        return temperature - self.atm.delta_temp

    @typechecked
    def f_ebar(self, temperature: xr.DataArray, rh_loc: xr.DataArray) -> xr.DataArray:
        """Flux e_bar.

        Args:
            temperature (xr.DataArray): temperature in kelvin.
            rh_loc (xr.DataArray): Relative humidity (dimemsionless)

        Returns:
            xr.DataArray: evaporation.
        """
        q_a = rh_loc * self.f_qs(temperature)

        # q_a is the surface-specific humidity
        return q_a * self.atm.p_s / self.atm.e_factor

    # ------------ heat flux functions:

    @typechecked
    def f_qlw1(
        self,
        temperature: xr.DataArray,
        cloud_cover: Union[xr.DataArray, float],
        f: float,
        rh_loc: xr.DataArray,
    ) -> xr.DataArray:
        """The first term of the long wave flux equation (14).

        Qlw1 = epsilon sigma T^4 f' (1 - a_cloud_const C^2)

        Args:
            temperature (xr.DataArray): temperature of the surface in kelvin.
            cloud_cover (Union[xr.DataArray, float]): the cloud cover, which could be
                a constant or an array.
            f (float):
            rh_loc (xr.DataArray): relative humidity

        Returns:
            xr.DataArray: The first term of the long wave flux equation.
        """
        temp_a = self.f_temp_a(temperature)
        if not self.atm.vary_cloud_const:
            a_cloud_const = self.atm.a_cloud_const
        else:
            a_cloud_const = self.get_cloud_const(temperature)

        return (
            self.atm.qlw_coeff
            * (1 - a_cloud_const * cloud_cover ** 2)
            # bar(Ts)^4
            * temp_a ** 4
            * (f - self.atm.f2 * np.sqrt(self.f_ebar(temperature, rh_loc)))
            # f1'
        )

    def get_cloud_const(self, temperature: xr.DataArray) -> xr.DataArray:
        """
        Function to produce variable cloud constant using temperature.

        Models an approximation to the effect of deep convection.

        Args:
            temperature (xr.DataArray): temperature in Kelvin.

        Returns:
            xr.DataArray: The a constant. Normally 0.4 or 0.8. A standard value of
                0.6 is applied if you disable deep convection.
        """
        temperature = temperature - zero_Celsius
        cloud_const_da = temperature.copy()
        cloud_const_da = cloud_const_da.where(
            temperature >= self.atm.dc_threshold_temp
        ).fillna(self.atm.a_cloud_const_dc)
        cloud_const_da = cloud_const_da.where(
            temperature < self.atm.dc_threshold_temp
        ).fillna(self.atm.a_cloud_const_norm)
        return cloud_const_da

    @typechecked
    def f_qlw2(self, temperature: xr.DataArray) -> xr.DataArray:
        """Second term in QLW equation (14).

        Args:
            temperature (xr.DataArray): temperature in Kelvin.

        Returns:
            xr.DataArray: heat from long wave heat flux term 2.

        """
        return (
            4
            * self.atm.emmisivity
            * self.atm.stefan_boltzman_const
            * temperature ** 3
            * (temperature - self.f_temp_a(temperature))
        )

    @typechecked
    def f_qlw(
        self,
        temp: xr.DataArray,
        cloud_cover: xr.DataArray,
        f: float,
        rh_loc: xr.DataArray,
    ) -> xr.DataArray:
        """Full long wave heat flux.

        Args:
            temp (xr.DataArray): temperature in kelvin.
            cloud_cover (xr.DataArray): cloud cover, float.
            f (float): f.
            rh_loc (xr.DataArray): relative humidity.

        Returns:
            xr.DataArray: Full long wave heat flux.
        """
        return self.f_qlw1(temp, cloud_cover, f, rh_loc) + self.f_qlw2(temp)

    @typechecked
    def f_dqlw_df(
        self, temperature: xr.DataArray, cloud_cover: Union[xr.DataArray, float]
    ) -> xr.DataArray:
        """Flux dqlw_df.

        Args:
            temperature (xr.DataArray): temperature in Kelvin.
            cloud_cover (Union[xr.DataArray, float]): constant.

        Returns:
            xr.DataArray: flux dqlw_df.
        """
        if not self.atm.vary_cloud_const:
            a_cloud_const = self.atm.a_cloud_const
        else:
            a_cloud_const = self.get_cloud_const(temperature)
        return (
            self.atm.qlw_coeff
            * (1 - a_cloud_const * cloud_cover ** 2)
            * temperature ** 4
        )

    @typechecked
    def f_dqlw_dtemp(
        self,
        temperature: xr.DataArray,
        cloud_cover: Union[xr.DataArray, float],
        f: float,
        rh_loc: xr.DataArray,
    ) -> xr.DataArray:
        """Flux dqlw_dtemp.

        Args:
            temperature (xr.DataArray): Temperature dataarray in Kelvin.
            cloud_cover (Union[xr.DataArray, float]): Cloud cover.
            f (float): [description]
            rh_loc (xr.DataArray): the relative humidity.

        Returns:
            xr.DataArray: [description]

        """
        e_bar = self.f_ebar(temperature, rh_loc)
        q_s = self.f_qs(temperature)
        # q_a is the surface-specific humidity
        # q_s(Ts) is the saturation-specific humidity at the SST
        dqs_dtemp = self.f_dqs_dtemp(temperature)
        if not self.atm.vary_cloud_const:
            a_cloud_const = self.atm.a_cloud_const
        else:
            a_cloud_const = self.get_cloud_const(temperature)
        return self.atm.qlw_coeff * (
            (1 - a_cloud_const * cloud_cover ** 2)
            * temperature ** 3
            * (
                4 * f
                - self.atm.f2 * np.sqrt(e_bar) * (4 + temperature * dqs_dtemp / 2 / q_s)
            )
            + 12 * temperature ** 2 * self.atm.delta_temp
        )

    @typechecked
    def f_qa(self, t_s: np.ndarray, s_p: np.ndarray) -> np.ndarray:
        """Flux qa.

        Args:
            t_s (np.ndarray): sst in Kelvin.
            s_p (np.ndarray): surface pressure in mb.

        Returns:
            np.ndarray: q_a, surface specific humidity.

        """
        e_s = self.atm.es_0 * np.exp(
            17.67 * (t_s - zero_Celsius) / ((t_s - zero_Celsius) + zero_Celsius)
        )
        return self.atm.e_factor * self.atm.relative_humidity * e_s / s_p

    @typechecked
    def f_qa2(self, temp_surface: np.ndarray) -> np.ndarray:
        """flux qa2.

        Args:
            temp_surface (np.ndarray): sst in Kelvin.

        Returns:
            np.ndarray: q_s, surface specific humidity.

        """
        return 0.001 * (temp_surface - zero_Celsius - 11.0)

    @typechecked
    def f_evap(self, mask: np.ndarray, q_a: np.ndarray, wnsp: np.ndarray) -> np.ndarray:
        """evaporation flux.

        Args:
            mask (np.ndarray): land mask
            q_a (np.ndarray): surface air humidity
            wnsp (np.ndarray): surface windspeed in m/s

        Returns:
            np.ndarray: Evap in kg/m^2/s.

        """
        c_s_e = 0.0015 * (1 + mask / 2)
        # c_s_e = 0.0012
        return (
            c_s_e
            * self.atm.rho_air
            * (1 - self.atm.relative_humidity)
            * q_a
            * wnsp
            / self.atm.relative_humidity
        )

    @typechecked
    def f_mc(self, q_a: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """moisture convergence flux.

        To calculate surface heat fluxes and atmospheric moisture
        convergence, relative humidity is assumed to be spatially
        uniform in our standard model.
        (N.B., v is on y_axis_v points, u,q are on y_axis_u points)

        Args:
            q_a (np.ndarray): surface air humidity
            u (np.ndarray): low level winds in m/s
            v (np.ndarray): low level winds in m/s

        Returns:
            np.ndarray: Moisture convergence in kg/m^2/s.

        """
        qu = q_a * u
        qux = ifft(1.0j * self.kk_wavenumber * fft(qu) / R_earth.to_value()).real
        aq = (q_a[1 : self.atm.ny - 1, :] + q_a[0 : self.atm.ny - 2, :]) / 2.0
        qv = aq * v[1 : self.atm.ny - 1, :]
        z = np.zeros((1, self.atm.nx))
        qv = np.concatenate((z, qv, z), axis=0)
        # qvy = qv.diff('Yu')/dym
        qvy = (qv[1 : self.atm.ny, :] - qv[0 : self.atm.ny - 1, :]) / self.dym
        return -self.atm.h_q * (qux + qvy) * self.atm.rho_air

    # ---------------- equation solvers ---------------------

    @typechecked
    def tdma_solver(
        self,
        ny_loc: int,
        a_loc: np.ndarray,
        b: np.ndarray,
        c: np.ndarray,
        d: np.ndarray,
    ) -> np.ndarray:
        """tdma solver.

        'tdma_solver'  0.00243 s

        Tri Diagonal Matrix Algorithm (a.k.a. Thomas algorithm) solver

        E.g.
        https://gist.github.com/cbellei/8ab3ab8551b8dfc8b081c518ccd9ada9

        Args:
            ny_loc (int): local version of ny_loc.
            a_loc (np.ndarray): [description]
            b (np.ndarray): [description]
            c (np.ndarray): [description]
            d (np.ndarray): [description]

        Returns:
            np.ndarray: xc

        """

        nf = ny_loc  # number of equations
        ac, bc, cc, dc = map(np.array, (a_loc, b, c, d))  # copy arrays

        for it in range(1, nf):
            mc = ac[it, :] / bc[it - 1, :]
            bc[it, :] = bc[it, :] - mc * cc[it - 1, :]
            dc[it, :] = dc[it, :] - mc * dc[it - 1, :]

        xc = bc
        xc[-1, :] = dc[-1, :] / bc[-1, :]

        for il in range(nf - 2, -1, -1):
            xc[il, :] = (dc[il, :] - cc[il, :] * xc[il + 1, :]) / bc[il, :]

        return xc

    @typechecked
    def s91_solver(self, q1: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """S91 solver from TCAM.py.

        Presumably name refers to Richard Seager 1991.

        's91_solver'  0.00564 s

        The atmosphere equations are solved by Fourier transforming in longitude,
        forming an equation for v for each zonal wavenumber that is finite
        differenced, and the resulting tri-diagonal system is solved by matrix
        inversion, transforming back into longitude.

        Usef fft, ifft.

              g . pi . N ^ 2
        q1 = ---------------- . (k . theta_s . Q_c)
              theta_00 . z_t

        Args:
            q1 (np.ndarray): modified heating that drives winds.

        Returns:
            Tuple[np.ndarray, np.ndarray, np.ndarray]: u, v, phi

        """

        q1_time = fft(q1)
        f_q = self.fcu[:, np.newaxis] * q1_time
        a_f_q = (f_q[1 : self.atm.ny - 1, :] + f_q[0 : self.atm.ny - 2, :]) / 2.0
        km = self.kk_wavenumber / R_earth.to_value()
        d_q = (
            q1_time[1 : self.atm.ny - 1, :] - q1_time[0 : self.atm.ny - 2, :]
        ) / self.dym
        rk = (
            1.0j * km * self.atm.beta
            - self.atm.eps_u * self.atm.eps_v * self.atm.eps_p
            - self.atm.eps_v * km ** 2
        )

        fcp = self.fcu[1 : self.atm.ny - 1] ** 2 / 4.0
        fcm = self.fcu[0 : self.atm.ny - 2] ** 2 / 4.0

        ak = self.atm.eps_u / self.dym_2 - self.atm.eps_p * fcm[:, np.newaxis]
        ck = self.atm.eps_u / self.dym_2 - self.atm.eps_p * fcp[:, np.newaxis]
        bk = (
            -2 * self.atm.eps_u / self.dym_2
            - self.atm.eps_p * (fcm[:, np.newaxis] + fcp[:, np.newaxis])
            + rk[np.newaxis, :]
        )
        dk = -self.atm.eps_u * d_q + 1.0j * km[np.newaxis, :] * a_f_q

        # find tdma using it
        vtk = self.tdma_solver(self.atm.ny - 2, ak, bk, ck, dk)

        z = np.zeros((1, self.atm.nx))
        v_t = np.concatenate((z, vtk, z), axis=0)
        av = (v_t[1 : self.atm.ny, :] + v_t[0 : self.atm.ny - 1, :]) / 2.0
        fav = self.fcu[:, np.newaxis] * av
        dv = (v_t[1 : self.atm.ny, :] - v_t[0 : self.atm.ny - 1, :]) / self.dym
        coeff = self.atm.eps_u * self.atm.eps_p + km * km
        u_t = (
            self.atm.eps_p * fav + 1.0j * (q1_time + dv) * km[np.newaxis, :]
        ) / coeff[np.newaxis, :]
        phi_t = -(q1_time + 1.0j * u_t * km[np.newaxis, :] + dv) / self.atm.eps_p
        v = ifft(v_t).real
        u = ifft(u_t).real
        phi = ifft(phi_t).real

        return (u, v, phi)

    # -------------- smoother --------------------------

    # pylint: disable=dangerous-default-value
    @typechecked
    def smooth121(
        self,
        da: xr.DataArray,
        sdims: list,
        number_smooths: int = 1,
        perdims: list = list(),
    ) -> xr.DataArray:
        """Applies [0.25, 0.5, 0.25] stencil in sdims, one at a time.

        Args:
            da (xr.DataArray): xarray.DataArray - e.g., ds.var
            sdims (list): list of dimensions over which to smooth - e.g., ['lat','lon']
            number_smooths (int, optional): integer number of
                smooths to apply - e.g., 1. Defaults to 1.
            perdims (list, optional): list of dimension to be treated as period
                boundaries - e.g., ['lon']. Defaults to [].

        Returns:
            xr.DataArray: smoothed output.

        """
        mask = da.notnull()
        weight = xr.DataArray([0.25, 0.5, 0.25], dims=["window"])
        v = da.copy()
        origdims = v.dims

        for dim in sdims:
            for _ in range(0, number_smooths):
                if dim in perdims:
                    v0 = xr.concat(
                        [v.isel(**{dim: -1}), v, v.isel(**{dim: 0})], dim=dim
                    )
                else:
                    v0 = xr.concat(
                        [v.isel(**{dim: 0}), v, v.isel(**{dim: -1})], dim=dim
                    )
                v1 = v0.bfill(dim, limit=1)
                v0 = v1.ffill(dim, limit=1)
                v1 = v0.rolling(**{dim: 3}, center=True).construct("window").dot(weight)
                v = v1.isel(**{dim: slice(1, -1, None)})

        return v.where(mask, np.nan).transpose(*origdims)

    # ------------------ output functions -------------------------

    @timeit
    @typechecked
    def output_trends(self) -> None:
        """Output trends ds.

        Runs the Matsuno-Gill model with the trends in preipitation and
        sea surface temperature half added and half subtracted to work out the
        change in the other variables.

        𝜀𝑢 . 𝑢 − 𝑓 . 𝑣 + 𝜙 . 𝑥 =  0   (1)

        𝜀𝑣 . 𝑣 + 𝑓 . 𝑢 + 𝜙 . 𝑦 =  0   (2)

        𝜀𝜙 . 𝜙 + 𝑢 . 𝑥 + 𝑣 . 𝑦 = −𝑄1  (3)

        """

        ds = xr.Dataset(
            {
                "X": ("X", self.x_axis),
                "Yu": ("Yu", self.y_axis_u),
                "Yv": ("Yv", self.y_axis_v),
            }
        )
        ds.X.attrs = [("units", "degree_east")]
        ds.Yu.attrs = [("units", "degree_north")]
        ds.Yv.attrs = [("units", "degree_north")]

        ds["K"] = self.atm.k_days
        ds.K.attrs = [("units", "day")]
        ds["epsu"] = self.atm.eps_days
        ds.epsu.attrs = [("units", "day")]
        ds["epsv"] = self.atm.eps_days / self.atm.e_frac
        ds.epsv.attrs = [("units", "day")]
        ds["hq"] = self.atm.h_q
        ds.hq.attrs = [("units", "m")]

        # CLIMATOLOGIES

        def get_clim():
            # the average condions from ECMWF
            # Gets the windspeed, surface temperature, precipation, and surface pressure
            # VAR_DICT = {0: "ts", 1: "clt", 2: "sfcWind", 3: "rh", 4: "pr", 5: "ps", 6: "tau"}
            # "sfcWind"
            ds_clim = xr.open_dataset(self.setup.clim_file("sfcWind", path=True))
            # this needs to be replaced
            fwnsp = interp2d(ds_clim.X, ds_clim.Y, ds_clim.sfcWind, kind="linear")
            ds_clim = xr.open_dataset(self.setup.ts_clim(self.it))
            fts = interp2d(ds_clim.X, ds_clim.Y, ds_clim.ts, kind="linear")
            ds_clim = xr.open_dataset(self.setup.clim_file("pr", path=True))
            fpr = interp2d(ds_clim.X, ds_clim.Y, ds_clim.pr, kind="linear")
            ds_clim = xr.open_dataset(self.setup.clim_file("ps", path=True))
            fps = interp2d(ds_clim.X, ds_clim.Y, ds_clim.ps, kind="linear")
            # Return interpolation objects
            return fwnsp, fts, fpr, fps

        fwnsp, fts, fpr, fps = get_clim()

        # Regrid:

        wnsp = fwnsp(self.x_axis, self.y_axis_u)
        wnsp[wnsp < self.atm.wnsp_min] = self.atm.wnsp_min
        ds["wnspClim"] = (["Yu", "X"], wnsp)
        ds["tsClim"] = (["Yu", "X"], fts(self.x_axis, self.y_axis_u))
        ds["prClim"] = (["Yu", "X"], fpr(self.x_axis, self.y_axis_u))
        ds["spClim"] = (["Yu", "X"], fps(self.x_axis, self.y_axis_u))

        # TRENDS
        def get_trend():
            # return interpolation functions of the original trends
            ds_trend = xr.open_dataset(self.setup.ts_trend(self.it))
            fts_trend = interp2d(ds_trend.X, ds_trend.Y, ds_trend.ts, kind="linear")
            ds_trend = xr.open_dataset(
                self.setup.clim_file("pr", typ="trend", path=True)
            )
            fpr_trend = interp2d(ds_trend.X, ds_trend.Y, ds_trend.pr, kind="linear")
            return fts_trend, fpr_trend

        fts_trend, fpr_trend = get_trend()

        ts_trend = fts_trend(self.x_axis, self.y_axis_u)
        ds["tsTrend"] = (["Yu", "X"], ts_trend)

        pr_trend = fpr_trend(self.x_axis, self.y_axis_u)
        pr_trend[abs(self.y_axis_u) > 25] = 0
        # get rid of anything above 25 degrees north/south
        pr_trend[pr_trend > 5e-5] = 5e-5
        ds["prTrend"] = (["Yu", "X"], pr_trend)
        ds["prTrend"] = self.smooth121(ds.prTrend, ["Yu", "X"], perdims=["X"])

        # save a picture of the trend
        do_plot = False
        if do_plot:
            ds.prTrend.plot()
            plt.savefig(os.path.join(self.setup.atmos_path, "prTrend.png"))
            plt.clf()

        dsmask = xr.open_dataset(
            os.path.join(self.setup.atmos_data_path, "mask-360x180.nc")
        )
        fmask = interp2d(dsmask.X, dsmask.Y, dsmask.mask, kind="linear")
        ds["mask"] = (["Yu", "X"], fmask(self.x_axis, self.y_axis_u))

        # tsClim = ds.tsClim.values
        sp_clim = ds.spClim.values
        wnsp_clim = ds.wnspClim.values
        wnsp_clim[wnsp_clim < self.atm.wnsp_min] = self.atm.wnsp_min
        mask = ds.mask.values
        w_end = wnsp_clim
        w_beg = wnsp_clim

        # ts trend only appplied where not masked?
        # this is the vital part where the differences are added/subracted
        ts_end = (ds.tsClim + (1 - mask) * ds.tsTrend / 2).values
        ts_beg = (ds.tsClim - (1 - mask) * ds.tsTrend / 2).values
        # only used for land?
        pr_end = (ds.prClim + ds.prTrend / 2).values
        pr_beg = (ds.prClim - ds.prTrend / 2).values

        # ts trend passed to q_th
        q_th_end = (
            self.atm.newtonian_cooling_coeff_k1 * (ts_end - 30) / self.atm.b_coeff
        )
        q_th_beg = (
            self.atm.newtonian_cooling_coeff_k1 * (ts_beg - 30) / self.atm.b_coeff
        )
        # passed to qa_end
        qa_end = self.f_qa(ts_end, sp_clim)
        # qa_end = f_qa2(ts_end)
        e_end = self.f_evap(mask, qa_end, wnsp_clim)
        pr_c_end = e_end
        pr_c_end[pr_c_end < 0] = 0
        # pr_end[pr_end>pr_max] = pr_max

        qa_beg = self.f_qa(ts_beg, sp_clim)
        # qa_beg = f_qa2(ts_beg)
        e_beg = self.f_evap(mask, qa_beg, wnsp_clim)
        # precipitation is initially the same as evaporation
        pr_c_beg = e_beg  # precipation had already been defined.
        pr_c_beg[pr_c_beg < 0] = 0
        # pr_beg[pr_beg>pr_max] = pr_max

        def iterate(pr, pr_c, q_th, e1, qa1) -> Tuple[np.ndarray]:
            """
            Iterate through.

            Args:
                pr ([type]): precipitation default.
                pr_c ([type]): new precipation.
                q_th ([type]): q th.
                e1 ([type]): evaporation.
                qa1 ([type]): heat flux.

            Returns:
                Tuple[np.ndarray]: pr_c, u1, v1, phi1, mc1
            """
            # Find total pr, u and v at end
            for _ in range(0, self.atm.number_iterations):
                # Start main calculation
                q_c = (
                    np.pi
                    * self.atm.latent_heat_vap
                    * pr_c
                    / (
                        2
                        * self.atm.cp_air
                        * self.atm.rho_00
                        * self.atm.height_tropopause
                    )
                )  # heating from precip
                # convective heating part, Qc
                q1 = self.atm.b_coeff * (q_c + q_th)
                # Q1 is a modified heating since the part
                # involving θ is on the left-hand side
                (u1, v1, phi1) = self.s91_solver(q1)
                d_amc = xr.DataArray(self.f_mc(qa1, u1, v1), dims=["Yu", "X"])
                mc1 = self.smooth121(d_amc, ["Yu", "X"], perdims=["X"]).values
                if self.atm.prcp_land:
                    pr_c = (1 - mask) * (mc1 + e1) + mask * pr
                else:
                    pr_c = (1 - mask) * (mc1 + e1)
                pr_c[pr_c < 0] = 0
                # pr[pr > pr_max] = pr_max
            return pr_c, u1, v1, phi1, mc1

        # pr, pr_c, q_th, e1, qa1
        # pr_c, u1, v1, phi1, mc1

        pr_c_end, u_end, v_end, phi_end, mc_end = iterate(
            pr_end,
            pr_c_end,
            q_th_end,
            e_end,
            qa_end,
        )

        pr_c_beg, u_beg, v_beg, phi_beg, mc_beg = iterate(
            pr_beg,
            pr_c_beg,
            q_th_beg,
            e_beg,
            qa_beg,
        )

        # save and plot the trends
        ds["utrend"] = (["Yu", "X"], u_end - u_beg)
        ds["vtrend"] = (["Yv", "X"], v_end - v_beg)
        ds["phitrend"] = (["Yu", "X"], phi_end - phi_beg)
        ds["tstrend"] = (["Yu", "X"], ts_end - ts_beg)
        ds["PRtrend"] = (["Yu", "X"], pr_c_end - pr_c_beg)
        ds["Qthtrend"] = (["Yu", "X"], q_th_end - q_th_beg)
        ds["uend"] = (["Yu", "X"], u_end)  # u at end
        ds["vend"] = (["Yv", "X"], v_end)  # v at end.
        ds["wend"] = (["Yu", "X"], w_end)
        ds["phiend"] = (["Yu", "X"], phi_end)
        ds["tsend"] = (["Yu", "X"], ts_end)
        ds["PRend"] = (["Yu", "X"], pr_c_end)
        ds["Qthend"] = (["Yu", "X"], q_th_end)
        ds["Eend"] = (["Yu", "X"], e_end)
        ds["MCend"] = (["Yu", "X"], mc_end)
        ds["qaend"] = (["Yu", "X"], qa_end)
        ds["ubeg"] = (["Yu", "X"], u_beg)
        ds["vbeg"] = (["Yv", "X"], v_beg)
        ds["wbeg"] = (["Yu", "X"], w_beg)
        ds["phibeg"] = (["Yu", "X"], phi_beg)
        ds["tsbeg"] = (["Yu", "X"], ts_beg)
        ds["PRbeg"] = (["Yu", "X"], pr_c_beg)
        ds["Qthbeg"] = (["Yu", "X"], q_th_beg)
        ds["Ebeg"] = (["Yu", "X"], e_beg)
        ds["MCbeg"] = (["Yu", "X"], mc_beg)
        ds["qabeg"] = (["Yu", "X"], qa_beg)

        # There is 2 gridpoint noise in the phi field - so add a smooth in X:
        ds["phitrend"] = self.smooth121(
            ds.phitrend, ["X"], number_smooths=1, perdims=["X"]
        )

        do_plot = False

        if do_plot:

            # next ds
            ds_subset = ds.sel(X=slice(120, 290), Yu=slice(-40, 40))
            plt.figure(figsize=(8, 5))
            plt.subplot(311)
            ds_subset.utrend.plot(cmap="RdBu_r")
            plt.subplot(312)
            ds_subset.vtrend.plot(cmap="RdBu_r")
            plt.subplot(313)
            ds_subset.PRtrend.plot(cmap="RdBu_r")
            plt.savefig(
                os.path.join(self.setup.atmos_path, "S90-H2000-Stab.eps"),
                format="eps",
                dpi=1000,
            )
            plt.clf()

        # adding units.

        # these seem to be the most important quantities.

        ds.utrend.attrs = [("units", "m/s")]
        ds.vtrend.attrs = [("units", "m/s")]
        ds.phitrend.attrs = [("units", "m2/s2")]
        ds.PRtrend.attrs = [("units", "m/s")]
        ds.Qthtrend.attrs = [("units", "K/s")]

        en_dict = {
            "K": {"dtype": "f4"},
            "epsu": {"dtype": "f4"},
            "epsv": {"dtype": "f4"},
            "hq": {"dtype": "f4"},
        }

        outfile = self.setup.tcam_output()

        print(outfile)

        ds.to_netcdf(outfile, encoding=en_dict)

        # warnings.filterwarnings("ignore")
        do_plot = False  # False
        # vtrend = ds.vtrend.values
        # utrend = ds.utrend.values
        if do_plot:
            ftitle = (
                r"Winds: (u,v) from sst: "
                + r"  $K_1=1/"
                + str(self.atm.k_days)
                + r"$, $\epsilon=1/"
                + str(self.atm.eps_days)
                + r"$"
            )
            pr_trend = ds.prTrend.values

            # nsy = 2  # plot every nsy grid point
            # nsx = nsy
            plt.figure(figsize=(8, 8))
            plt.suptitle(ftitle, size=14)
            plt.subplot(211)
            plt.title(r"$(u,v)$ vectors, magnitude ($m/s$) in contours")
            # m.fillcontinents(color="grey")
            # Av = np.squeeze(vtrend[1:ny, :] + vtrend[0 : ny - 1, :]) / 2.0
            # Au = np.squeeze(utrend[:, :])
            # AQ = np.squeeze(pr_trend[:, :])
            # mag = np.sqrt(Au * Au + Av * Av)
            # CS = plt.contour(X, Yu, mag, 15, colors="k",
            #                  linewidths=0.2, vmin=0, vmax=10)
            # plt.clabel(CS, inline=1, fontsize=10, fmt="%.2f")
            # CS = plt.contour(X2,Yu,mag,5,linewidths=1)
            # plt.clabel(CS, inline=1, fontsize=10,fmt='%.1f')
            # print(np.shape(Au[::nsy,::nsx]),np.shape(Av[::nsy,::nsx]))
            print(ds)
            print(ds.utrend)
            print(ds.vtrend)
            ds.plot.quiver(x="X", y="Yu", u="utrend", v="vtrend")
            # m.quiver(X[::nsx], Yu[::nsy], Au[::nsy, ::nsx]
            # Av[::nsy, ::nsx], scale=75)
            plt.subplot(212)
            plt.title(r"$Qc (mm/day)$")
            # m.fillcontinents(color="grey")
            # m.pcolormesh(X, Yu, 24 * 3600 * AQ, cmap="RdBu_r", vmin=-5, vmax=5)
            # CS = plt.contour(
            #    x_axis, y_axis_u, 24 * 3600 * AQ, 15, colors="k",
            #  # linewidths=0.2, vmin=-20, vmax=20
            # )
            # plt.clabel(CS, inline=1, fontsize=10, fmt="%.2f")
            plt.tight_layout()
            plt.subplots_adjust(top=0.90)
            plt.savefig(
                os.path.join(
                    self.setup.atmos_path,
                    "windsFromSST-K"
                    + str(self.atm.k_days)
                    + "-eps"
                    + str(self.atm.eps_days)
                    + ".eps",
                ),
                format="eps",
                dpi=1000,
            )
            plt.clf()

    # ##--------------------------- Begin dQ ----------------------------

    def load_clim60(self) -> xr.Dataset:
        """
        Load the inputs to get_dclim before processing.

        This is ok for the first iteration, but will need to be changed.

        Returns:
            xr.Dataset: An mfdataset with "ts", "clt", "sfcWind", "rh".
        """

        files = []

        for i, m in enumerate(self.atm.mem):
            if i < 4:
                name = self.names[m]
                variable = self.var[i]
                if variable == "ts":
                    # the surface temperature can be an input from the ocean model.
                    file = self.setup.ts_clim60(self.it)
                    # temperature is in degrees kelvin!
                else:
                    file = self.setup.clim60_name(i, path=True)
                assert os.path.isfile(file)
                print(name, variable, file)
                files += [file]  # append to list.

        return xr.open_mfdataset(files, decode_times=False)

    @timeit
    @typechecked
    def get_dclim(
        self,
    ) -> Tuple[
        xr.Dataset,
        xr.DataArray,
        xr.DataArray,
        xr.DataArray,
        xr.DataArray,
        xr.DataArray,
        xr.DataArray,
        xr.DataArray,
        xr.DataArray,
    ]:
        """Opens the files, and applies functions to get surface fluxes.

        Get the surface fluxes, qd_df, dq_dT among other things.

        Example:
            Logic behind calculating the flux parameters::

                dq_dt = self.f_dqlh_dtemp(t_sb, u_b, rh_b)
                        + self.f_dqlw_dtemp(t_sb, c_b, self.atm.f1_bar, rh)

                dq_df = self.f_dqlw_df(t_sb, c_b)

                dq_dt = (${atm.rho_air} * ${atm.c_e} * ${atm.latent_heat_vap}* u_sp
                            * (self.atm.e_factor * self.atm.es_0
                            * np.exp(17.67
                                * (temperature - self.atm.temp_0_c)
                                / (temperature - self.atm.temp_0_c + 243.5)
                            ) / self.atm.p_s
                            * (17.67 * 243.5)
                            / (temperature - self.atm.temp_0_c + 243.5) ** 2
                        ) * (1 - rh)
                        + ${atm.emmisivity} * ${atm.stefan_boltzman_const} * (
                            (1 - a_cloud_const * cloud_cover ** 2)
                            * temperature ** 3
                            * (
                                4 * f
                                - self.atm.f2 * np.sqrt(e_bar)
                                * (4 + temperature * dqs_dtemp / 2 / q_s)
                            )
                            + 12 * temperature ** 2 * self.atm.delta_temp
                            )
                        )

                dq_df = (
                            ${atm.emmisivity} * ${atm.stefan_boltzman_const}
                            * (1 - a_cloud_const * cloud_cover ** 2)
                            * temperature ** 4
                        )


        Returns:
            Tuple[
                xr.Dataset,
                xr.DataArray,
                xr.DataArray,
                xr.DataArray,
                xr.DataArray,
                xr.DataArray,
                xr.DataArray,
                xr.DataArray,
                xr.DataArray,
            ]: dclim, u_b, alh, alw, blw, dtemp_se, rh, c_b, t_sb.

        """
        # set Q'_LW + Q'_LH = 0, solve for Ts' (assuming U'=0)
        # Q'_LW = (alw(temp_surface_bar,c_bar,f1_bar)* Tsprime
        #          + blw(temp_surface_bar,c_bar) * f1prime)
        # Q'_LH = alh(temp_surface_bar,u_bar) * Tsprime
        # Q'_LH is from formula 13 in paper
        # Q'_LW is from formula 14 in paper

        dclim_loc = self.load_clim60()
        t_sb_loc = 1.0 * dclim_loc.ts
        # process the climatological windspeed
        tmp_wsp = 1.0 * dclim_loc.sfcWind.stack(z=("lon", "lat")).load()
        # clip windspeed to above 4 ms-1
        tmp_wsp[tmp_wsp < self.atm.wnsp_min] = self.atm.wnsp_min
        u_b_loc = tmp_wsp.unstack("z").T  # climatological windspeed.
        c_b_loc = dclim_loc.clt / 100.0
        rh_loc = dclim_loc.rh / 100.0
        f1p = self.atm.f1prime  # -0.003  # f1prime

        alh0 = self.f_dqlh_dtemp(t_sb_loc, self.atm.u_bar, rh_loc)
        alw0 = self.f_dqlw_dtemp(t_sb_loc, self.atm.c_bar, self.atm.f1_bar, rh_loc)
        blw0 = self.f_dqlw_df(t_sb_loc, self.atm.c_bar)
        dtemp_se0 = -blw0 * f1p / (alh0 + alw0)
        dclim_loc["dTse0"] = dtemp_se0

        alh1 = self.f_dqlh_dtemp(t_sb_loc, u_b_loc, rh_loc)
        alw1 = self.f_dqlw_dtemp(t_sb_loc, self.atm.c_bar, self.atm.f1_bar, rh_loc)
        blw1 = self.f_dqlw_df(t_sb_loc, self.atm.c_bar)
        dtemp_se1 = -blw1 * f1p / (alh1 + alw1)
        dclim_loc["dTse1"] = dtemp_se1

        alh2 = self.f_dqlh_dtemp(t_sb_loc, self.atm.u_bar, rh_loc)
        alw2 = self.f_dqlw_dtemp(t_sb_loc, c_b_loc, self.atm.f1_bar, rh_loc)
        blw2 = self.f_dqlw_df(t_sb_loc, c_b_loc)
        dtemp_se2 = -blw2 * f1p / (alh2 + alw2)
        dclim_loc["dTse2"] = dtemp_se2

        alh_loc = self.f_dqlh_dtemp(t_sb_loc, u_b_loc, rh_loc)
        alw_loc = self.f_dqlw_dtemp(t_sb_loc, c_b_loc, self.atm.f1_bar, rh_loc)
        blw_loc = self.f_dqlw_df(t_sb_loc, c_b_loc)
        dtemp_se_loc = -blw_loc * f1p / (alh_loc + alw_loc)

        dclim_loc["dTse"] = dtemp_se_loc
        dclim_loc["ALH"] = alh_loc
        dclim_loc["ALW"] = alw_loc
        dclim_loc["BLW"] = blw_loc
        dclim_loc["QLW"] = alw_loc + blw_loc * f1p / dtemp_se_loc

        dclim_loc.to_netcdf(self.setup.q_output())

        return (
            dclim_loc,
            u_b_loc,
            alh_loc,
            alw_loc,
            blw_loc,
            dtemp_se_loc,
            rh_loc,
            c_b_loc,
            t_sb_loc,
        )

    @typechecked
    def output_dq(self) -> None:
        """Outputs "dQ.nc"."""

        dclim, u_b, alh, alw, blw, dtemp_se, rh, c_b, t_sb = self.get_dclim()

        # Now, save the dq_df and dq_dt terms for using in TCOM:
        dq_dt = alh + alw
        dq_df = blw

        # Define the new Dataset
        dq = xr.Dataset(
            {
                "lon": ("lon", dclim["lon"].values),
                "lat": ("lat", dclim["lat"].values),
                "dq_dt": (["lat", "lon"], dq_dt.data),
                "dq_df": (["lat", "lon"], dq_df.data),
            }
        )

        dq.lon.attrs = dclim.lon.attrs
        dq.lat.attrs = dclim.lat.attrs
        dq.dq_dt.attrs = [("units", "W/m^2/K")]
        dq.dq_df.attrs = [("units", "W/m^2")]

        dq["ALH"] = alh
        dq["ALW"] = alw
        dq["BLW"] = blw
        dq["dTse"] = dtemp_se
        dq["rh"] = rh
        dq["Ub"] = u_b
        dq["Cb"] = c_b
        dq["Tsb"] = t_sb

        dq.to_netcdf(self.setup.dq_output())

    @typechecked
    def make_figure(
        self,
        cmap: Union[str, Any] = "viridis",
        lat: str = "latitude",
        lon: str = "longitude",
    ) -> None:
        """Make figure.

        Args:
            cmap (str, optional): matplotlib colormap. Defaults to "viridis".
            lat (str, optional): latitude label name. Defaults to "latitude".
            lon (str, optional): longitude label name. Defaults to "longitude".

        """

        dclim, u_b, _, _, _, _, _, _, _ = self.get_dclim()

        v_min_ad = 0.0
        v_max_ad = 0.6

        plt.figure(figsize=(8, 6))
        plt.subplot(321)
        dp = dclim.dTse0.plot.contourf(
            levels=11, cmap=cmap, vmin=v_min_ad, vmax=v_max_ad, add_colorbar=0
        )
        plt.title(r"$T^{\,\prime}_s$  for $\bar U,\bar C$")
        plt.ylabel(lat)
        plt.xlabel(lon)
        _ = plt.colorbar(dp)
        plt.subplot(322)
        dp = dclim.dTse1.plot.contourf(
            levels=11, cmap=cmap, vmin=v_min_ad, vmax=v_max_ad, add_colorbar=0
        )
        plt.title(r"$T^{\,\prime}_s$  for $\bar U(x,y), \bar C$")
        plt.ylabel(lat)
        plt.xlabel(lon)
        _ = plt.colorbar(dp)
        plt.subplot(323)
        dp = dclim.dTse2.plot.contourf(
            levels=11, cmap=cmap, vmin=v_min_ad, vmax=v_max_ad, add_colorbar=0
        )
        plt.title(r"$T^{\,\prime}_s$  for $\bar U, \bar C(x,y)$")
        plt.ylabel(lat)
        plt.xlabel(lon)
        _ = plt.colorbar(dp)
        plt.subplot(324)
        dp = dclim.dTse.plot.contourf(
            levels=11, cmap=cmap, vmin=v_min_ad, vmax=v_max_ad, add_colorbar=0
        )
        plt.title(r"$T^{\,\prime}_s$  for $\bar U(x,y), \bar C(x,y)$")
        plt.ylabel(lat)
        plt.xlabel(lon)
        _ = plt.colorbar(dp)
        plt.subplot(325)
        dp = (dclim.clt / 100).plot.contourf(
            levels=11, cmap=cmap, vmin=0.0, vmax=1.0, add_colorbar=0
        )
        plt.title(r"$\bar C(x,y)$")
        plt.ylabel(lat)
        plt.xlabel(lon)
        _ = plt.colorbar(dp)
        plt.subplot(326)
        dp = u_b.plot.contourf(levels=11, cmap=cmap, vmin=4.0, vmax=8.0, add_colorbar=0)
        plt.title(r"$\bar U(x,y)$")
        plt.ylabel(lat)
        plt.xlabel(lon)
        _ = plt.colorbar(dp)
        plt.tight_layout()
        plt.savefig(
            os.path.join(self.setup.atmos_path, "Tsp4.eps"), format="eps", dpi=1000
        )

    def run_all(self, it: int = 0) -> None:
        self.it = it
        self.output_trends()
        self.output_dq()
        self.make_figure()
