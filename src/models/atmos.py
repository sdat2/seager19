"""src.models.atmos.

pytest src/test/test_atmos.py
python3 src/models/atmos.py

Model solution method:

The atmosphere equations are solved by Fourier transforming in longitude,
forming an equation for v for each zonal wavenumber that is finite
differenced, and the resulting tri-diagonal system is solved by matrix
inversion, transforming back into longitude. Finally, u and Φ are derived
by back-substitution. The ocean equations are solved using the ‘INC’
scheme31, integrating the model forward, after spin-up with
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

"""
from typing import Tuple, Union
import os
import numpy as np
from scipy.interpolate import interp2d
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt
import logging
import xarray as xr
from src.constants import ATMOS_TMP_PATH, ATMOS_DATA_PATH, ATMOS_PATH
from src.utils import timeit


log = logging.getLogger(__name__)


# ------------- constants -----------------------
# begining TCAM
k_days = 10  # K = 1/10
e_frac = 2.0  # multiply eps_u by e_frac to get eps_v
h_q = 1800  # Hq is a scale depth for moisture
prcp_land = 1  # use data precip trends over land
wnsp_min = 4
rho_00: float = 0.3
# this means the density at the surface.
# this is actual rho_bar in the paper
# 0.3 kg m-3
pr_max: float = 20.0 / 3600 / 24  # 20 / seconds in hour / hours in day.
relative_humidity: float = 0.80  # relative humidity uniformly 0.8
number_iterations: int = 50  #  int
gravity = 9.8  #  m s-2
height_tropopause = 15e3  # metres
theta_00 = 300  # potential temperature at the surface.
nbsq = 3.0e-4  # N^2 s-2. N^2 is a specified buoyancy frequency.
radius_earth = 6.37e6  # metres
sec_in_day = 86400  # seconds in day.
omega2 = 2 * (2 * np.pi / sec_in_day)  # 2 * rad per second
latent_heat_vap = 2.5e6  # latent heat # J kg-1
cp_air = 1000  #  cp_air is the specific heat capacity of air.
# J kg-1 K-1
b_coeff = gravity * np.pi / (nbsq * theta_00 * height_tropopause)
eps_days = 0.75
# Over 1958–2017, the CO2 changed from ~300 to ~400 ppm, which would be about 0.75 W m−2
# eps days might be the efficiency of entrainment.
eps = 1.0 / (eps_days * sec_in_day)  # 1/.75 d
eps_u = eps  # 1/.75 d
eps_v = e_frac * eps  # e_frac=1/2 in paper
# Newtonian cooling, K
K1 = b_coeff / (k_days * sec_in_day)
epsp = (np.pi / height_tropopause) ** 2 / (nbsq * k_days * sec_in_day)
beta = omega2 / radius_earth
rho_air = 1.225  # kg m-3 - also called rho_00
c_e = 0.00125  # 1.25e-3 # cE is an exchange coefficient
emmisivity = 0.97  # problem here with second definintion.
stefan_boltzman_const = 5.67e-8
p_s = 1000  # pressure at the surface?
es_0 = 6.11
delta_temp = 1.0  #  ΔT = 1 K
f2 = 0.05  #  f2 = 0.05
# 'a' should decrease when deep convection happens above 28 degC
#  a = Ts-temp_0_c;a[a>28] = 40;a[a<=28] = 80;a = 0.01*a
a = 0.6  # this isn't the option used in the paper.

# basic parameters
temp_0_c = 273.15  # zero degrees
f1_bar = 0.39  #  f1 = 0.39
# f'1  is the anomaly in f1—a parameter that can be adjusted
# to control the variation in surface longwave radiation due
# to a change in CO2
u_bar = 5.0
temp_surface_bar = temp_0_c + 25  # 25C in Kelvin
c_bar = 0.6

# grid characteristics

nx: int = 180  # number of x grid boxes
ny: int = 60  # this seems like half the grid space
y_north_lim = 60  # upper lat limit
y_south_lim = -y_north_lim

# make grids
dx: float = 360 / nx
dy = (y_north_lim - y_south_lim) / ny
x_axis = np.linspace(0, 360 - dx, nx)
y_axis_v = np.linspace(y_south_lim + dy / 2, y_north_lim - dy / 2, ny)
y_axis_u = np.linspace(y_south_lim + dy, y_north_lim - dy, ny - 1)
y_axis_i = np.linspace(y_south_lim + 3 * dy / 2, y_north_lim - 3 * dy / 2, ny - 2)
x_spacing = x_axis[1] - x_axis[0]  # degrees
y_spacing = y_axis_v[1] - y_axis_v[0]  # degrees
dxm = x_spacing * radius_earth * np.pi / 180
dym = y_spacing * radius_earth * np.pi / 180
dym_2 = dym * dym

# need to have the correct ordering of the wave numbers for fft
num = nx
if num % 2 == 0:
    kk_wavenumber = np.asarray(
        list(range(0, num // 2)) + [0] + list(range(-num // 2 + 1, 0)), np.float64
    )
else:
    kk_wavenumber = np.asarray(
        list(range(0, (num - 1) // 2)) + [0] + list(range(-(num - 1) // 2, 0)),
        np.float64,
    )

# Find linearization of Q_LH (latent heating)
const1: float = rho_air * c_e * latent_heat_vap

mem: str = "EEEf"

# the different model names in a dict?
names: dict = {
    "E": "ECMWF",
    "F": "ECMWF-orig",
    "B": "CMIP5-39m",
    "C": "CMIP5",
    "D": "CMIP5-orig",
    "H": "HadGEM2",
    "f": "fixed",
    "e": "fixed78",
    "g": "fixed82",
    "W": "WHOI",
    "M": "MERRA",
    "I": "ISCCP",
}

var: dict = {0: "ts", 1: "clt", 2: "sfcWind", 3: "rh"}


# --------------- flux functions ----------------------


def f_cor(y: np.ndarray) -> np.ndarray:
    """Corriolis force coeff.

    omega2 = 2 * (2 * np.pi / sec_in_day) # 2 * rad per second

    # TODO is this an extra factor of 2?

    Args:
        y (np.ndarray): latitude

    Returns:
        np.ndarray: Corriolis force coeff.
    """
    return omega2 * y * np.pi / 180


fcu = f_cor(y_axis_u)


def f_es(temperature: xr.DataArray) -> xr.DataArray:
    """Flux es.

    Args:
        temperature (xr.DataArray): temp.

    Returns:
        xr.DataArray: Flux es.
    """
    return es_0 * np.exp(
        17.67 * (temperature - temp_0_c) / (temperature - temp_0_c + 243.5)
    )


def f_qs(temperature: xr.DataArray) -> xr.DataArray:
    """flux q_s.

    q_s(Ts) is the saturation-specific humidity at the SST
    rqs(Ts) = q_s(Ts)

    Args:
        temperature (xr.DataArray): temp.

    Returns:
        xr.DataArray: Flux q_s.
    """
    return 0.622 * f_es(temperature) / p_s


def f_dqs_dtemp(temperature: xr.DataArray) -> xr.DataArray:
    """flux dqs.

    Args:
        temperature (xr.DataArray): temp.

    Returns:
        xr.DataArray: flux q_s.
    """
    return f_qs(temperature) * (17.67 * 243.5) / (temperature - temp_0_c + 243.5) ** 2


def f_qlh(
    temperature: xr.DataArray, u_sp: xr.DataArray, rh_loc: xr.DataArray
) -> xr.DataArray:
    """heat flux from latent heat.

    It is assumed that the surface heat flux anomaly is
    dominated by longwave and latent heat fluxes and
    that the solar radiation does not change and the
    sensible heat flux anomaly is small.

    Args:
        temperature (xr.DataArray): [description]
        u_sp (xr.DataArray): [description]
        rh_loc (xr.DataArray): relative humidity

    Returns:
        xr.DataArray: flux qlh.

    """
    # print("u_sp", type(u_sp))
    return const1 * u_sp * f_qs(temperature) * (1 - rh_loc)


def f_dqlh_dtemp(
    temperature: xr.DataArray, u_sp: xr.DataArray, rh_loc: xr.DataArray
) -> xr.DataArray:
    """flux dqlh_dtemp.

    Args:
        temperature (xr.DataArray): temperature (in kelvin?).
        u_sp (xr.DataArray): u speed.
        rh_loc (xr.DataArray): relative humidity.

    Returns:
        xr.DataArray: flux dqlh_dtemp.
    """
    # print("u_sp", type(u_sp))
    return const1 * u_sp * f_dqs_dtemp(temperature) * (1 - rh_loc)


# Find linearization of Q_LW (longwave)
const2 = emmisivity * stefan_boltzman_const


def f_temp_a(temperature: xr.DataArray) -> xr.DataArray:
    """temperature anomaly.

    Delta temp = 1 in paper.

    Args:
        temperature (xr.DataArray): temp.

    Returns:
        xr.DataArray: temperature anomaly.
    """
    return temperature - delta_temp


def f_ebar(temperature: xr.DataArray, rh_loc: xr.DataArray) -> xr.DataArray:
    """flux e_bar.

    Args:
        temperature (xr.DataArray): [description]
        rh_loc (xr.DataArray): [description]

    Returns:
        xr.DataArray: [description]
    """
    q_a = rh_loc * f_qs(temperature)

    # q_a is the surface-specific humidity
    return q_a * p_s / 0.622


def f_qlw1(
    temperature: xr.DataArray, cloud_cover: xr.DataArray, f: float, rh_loc: xr.DataArray
) -> xr.DataArray:
    """The first term of the long wave flux equation (14).

    Qlw1 = epsilon sigma T^4 f' (1 - a C^2)

    Args:
        temperature (xr.DataArray): temperature of the surface?
        cloud_cover (xr.DataArray): [description]
        f (float): [description]
        rh_loc (xr.DataArray): [description]

    Returns:
        xr.DataArray: The first term of the long wave flux equation.
    """
    # print("rh_loc", type(rh_loc))
    temp_a = f_temp_a(temperature)
    return (
        const2
        * (1 - a * cloud_cover ** 2)
        # bar(Ts)^4
        * temp_a ** 4
        * (f - f2 * np.sqrt(f_ebar(temperature, rh_loc)))
        # f1'
    )


def f_qlw2(temperature: xr.DataArray) -> xr.DataArray:
    """Second term in QLW equation (14).

    Args:
        temperature (xr.DataArray): [description]

    Returns:
        xr.DataArray: [description]

    """
    # print("temperature", type(temperature))
    return (
        4
        * emmisivity
        * stefan_boltzman_const
        * temperature ** 3
        * (temperature - f_temp_a(temperature))
    )


# def f_QLW(T, f, rh_loc):
#     return f_qlw1(T, f, rh_loc) + f_qlw2(T)
# This function seemed to call a function in an incorrect way.


def f_dqlw_df(temperature: xr.DataArray, cloud_cover: xr.DataArray) -> xr.DataArray:
    """flux dqlw_df.

    Args:
        temperature (xr.DataArray): temp.
        cloud_cover (xr.DataArray): constant.

    Returns:
        xr.DataArray: flux dqlw_df.
    """
    return const2 * (1 - a * cloud_cover ** 2) * temperature ** 4


def f_dqlw_dtemp(
    temperature: xr.DataArray,
    cloud_cover: xr.DataArray,
    f: float,
    rh_loc: xr.DataArray
) -> xr.DataArray:
    """flux dqlw_dtemp.

    Args:
        temperature (xr.DataArray): [description]
        cloud_cover (xr.DataArray): [description]
        f (float): [description]
        rh_loc (xr.DataArray): [description]

    Returns:
        xr.DataArray: [description]

    """
    e_bar = f_ebar(temperature, rh_loc)
    q_s = f_qs(temperature)
    # q_a is the surface-specific humidity
    # q_s(Ts) is the saturation-specific humidity at the SST
    dqs_dtemp = f_dqs_dtemp(temperature)
    return const2 * (
        (1 - a * cloud_cover ** 2)
        * temperature ** 3
        * (4 * f - f2 * np.sqrt(e_bar)
        * (4 + temperature * dqs_dtemp / 2 / q_s))
        + 12 * temperature ** 2 * delta_temp
    )


def f_qa(ts: np.ndarray, sp: np.ndarray) -> np.ndarray:
    """f_qa.

    Args:
        ts (np.ndarray): sst in Kelvin
        sp (np.ndarray): surface pressure in mb

    Returns:
        np.ndarray: q_s, surface specific humidity
    """
    efac = 0.622
    es = es_0 * np.exp(17.67 * (ts - 273.15) / ((ts - 273.15) + 243.5))
    return efac * relative_humidity * es / sp


def f_qa2(temp_surface: np.ndarray) -> np.ndarray:
    """flux qa2.

    Args:
        temp_surface (np.ndarray): sst in Kelvin.

    Returns:
        np.ndarray: q_s, surface specific humidity.

    """
    return 0.001 * (temp_surface - 273.15 - 11.0)


def f_evap(mask: np.ndarray, q_a: np.ndarray, wnsp: np.ndarray) -> np.ndarray:
    """evaporation flux.

    Args:
        mask (np.ndarray): [description]
        q_a (np.ndarray): surface air humidity
        wnsp (np.ndarray): surface windspeed in m/s

    Returns:
        np.ndarray: Evap in kg/m^2/s.

    """
    c_s_e = 0.0015 * (1 + mask / 2)
    # c_s_e = 0.0012
    return c_s_e * rho_air * (1 - relative_humidity) * q_a * wnsp / relative_humidity


def f_mc(q_a: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
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
    qux = ifft(1.0j * kk_wavenumber * fft(qu) / radius_earth).real
    aq = (q_a[1 : ny - 1, :] + q_a[0 : ny - 2, :]) / 2.0
    qv = aq * v[1 : ny - 1, :]
    z = np.zeros((1, nx))
    qv = np.concatenate((z, qv, z), axis=0)
    # qvy = qv.diff('Yu')/dym
    qvy = (qv[1:ny, :] - qv[0 : ny - 1, :]) / dym
    return -h_q * (qux + qvy) * rho_air


# ---------------- equation solvers ---------------------


@timeit
def tdma_solver(
    ny_loc: int, a_loc: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """tdma solver.

    'tdma_solver'  0.00243 s

    Args:
        ny_loc (int):
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


@timeit
def s91_solver(q1: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """S91 folder from TCAM.py.

    's91_solver'  0.00564 s

     The atmosphere equations are solved by Fourier transforming in longitude,
     forming an equation for v for each zonal wavenumber that is finite
     differenced, and the resulting tri-diagonal system is solved by matrix
     inversion, transforming back into longitude.

     Usef fft, ifft.

               g . pi . N ^ 2
         q1 = ---------------- . (k theta_s Q_c)
               theta_00 . z_t

     Args:
         q1 (np.ndarray): modified heating that drives winds.

     Returns:
         Tuple[np.ndarray, np.ndarray, np.ndarray]: u, v, phi

    """
    # print("q1", type(q1))
    q1_time = fft(q1)
    f_q = fcu[:, np.newaxis] * q1_time
    a_f_q = (f_q[1 : ny - 1, :] + f_q[0 : ny - 2, :]) / 2.0
    km = kk_wavenumber / radius_earth
    d_q = (q1_time[1 : ny - 1, :] - q1_time[0 : ny - 2, :]) / dym
    rk = 1.0j * km * beta - eps_u * eps_v * epsp - eps_v * km ** 2
    fcp = fcu[1 : ny - 1] ** 2 / 4.0
    fcm = fcu[0 : ny - 2] ** 2 / 4.0

    ak = eps_u / dym_2 - epsp * fcm[:, np.newaxis]
    ck = eps_u / dym_2 - epsp * fcp[:, np.newaxis]
    bk = (
        -2 * eps_u / dym_2
        - epsp * (fcm[:, np.newaxis] + fcp[:, np.newaxis])
        + rk[np.newaxis, :]
    )
    dk = -eps_u * d_q + 1.0j * km[np.newaxis, :] * a_f_q

    vtk = tdma_solver(ny - 2, ak, bk, ck, dk)

    z = np.zeros((1, nx))
    vt = np.concatenate((z, vtk, z), axis=0)
    av = (vt[1:ny, :] + vt[0 : ny - 1, :]) / 2.0
    fav = fcu[:, np.newaxis] * av
    dv = (vt[1:ny, :] - vt[0 : ny - 1, :]) / dym
    coeff = eps_u * epsp + km * km
    ut = (epsp * fav + 1.0j * (q1_time + dv) * km[np.newaxis, :]) / coeff[np.newaxis, :]
    phit = -(q1_time + 1.0j * ut * km[np.newaxis, :] + dv) / epsp
    v = ifft(vt).real
    u = ifft(ut).real
    phi = ifft(phit).real
    return (u, v, phi)


# -------------- smoother --------------------------


# pylint: disable=dangerous-default-value
def smooth121(
    da: xr.DataArray,
    sdims: list,
    number_smooths: int = 1,
    perdims: list = list(),
) -> xr.DataArray:
    """Applies [0.25, 0.5, 0.25] stencil in sdims, one at a time.

    Args:
        da (xr.DataArray): xarray.DataArray - e.g., ds.var
        sdims (list): list of dimensions over which to smooth - e.g., ['lat','lon']
        number_smooths (int, optional): integer number of smooths to apply - e.g., 1.
            Defaults to 1.
        perdims (list, optional): list of dimension to be treated as period
            boundaries - e.g., ['lon']. Defaults to [].

    """
    mask = da.notnull()
    weight = xr.DataArray([0.25, 0.5, 0.25], dims=["window"])
    v = da.copy()
    origdims = v.dims

    for dim in sdims:
        for _ in range(0, number_smooths):
            if dim in perdims:
                v0 = xr.concat([v.isel(**{dim: -1}), v, v.isel(**{dim: 0})], dim=dim)
            else:
                v0 = xr.concat([v.isel(**{dim: 0}), v, v.isel(**{dim: -1})], dim=dim)
            v1 = v0.bfill(dim, limit=1)
            v0 = v1.ffill(dim, limit=1)
            v1 = v0.rolling(**{dim: 3}, center=True).construct("window").dot(weight)
            v = v1.isel(**{dim: slice(1, -1, None)})

    return v.where(mask, np.nan).transpose(*origdims)


# ------------------ output functions -------------------------


@timeit
def output_trends(direc: str = "") -> None:
    """output trends ds.

    𝜀𝑢 . 𝑢 − 𝑓 . 𝑣 + 𝜙 . 𝑥 =  0   (1)
    𝜀𝑣 . 𝑣 + 𝑓 . 𝑢 + 𝜙 . 𝑦 =  0   (2)
    𝜀𝜙 . 𝜙 + 𝑢 . 𝑥 + 𝑣 . 𝑦 = −𝑄1  (3)

    Args:
        direc (str): directory to save to.

    """
    log("Output trends.")

    ds = xr.Dataset(
        {"X": ("X", x_axis), "Yu": ("Yu", y_axis_u), "Yv": ("Yv", y_axis_v)}
    )
    ds.X.attrs = [("units", "degree_east")]
    ds.Yu.attrs = [("units", "degree_north")]
    ds.Yv.attrs = [("units", "degree_north")]

    ds["K"] = k_days
    ds.K.attrs = [("units", "day")]
    ds["epsu"] = eps_days
    ds.epsu.attrs = [("units", "day")]
    ds["epsv"] = eps_days / e_frac
    ds.epsv.attrs = [("units", "day")]
    ds["hq"] = h_q
    ds.hq.attrs = [("units", "m")]

    # CLIMATOLOGIES

    ds_clim = xr.open_dataset(os.path.join(ATMOS_DATA_PATH, "sfcWind-ECMWF-clim.nc"))
    fwnsp = interp2d(ds_clim.X, ds_clim.Y, ds_clim.sfcWind, kind="linear")
    ds_clim = xr.open_dataset(os.path.join(ATMOS_DATA_PATH, "ts-ECMWF-clim.nc"))
    fts = interp2d(ds_clim.X, ds_clim.Y, ds_clim.ts, kind="linear")
    ds_clim = xr.open_dataset(os.path.join(ATMOS_DATA_PATH, "pr-ECMWF-clim.nc"))
    fpr = interp2d(ds_clim.X, ds_clim.Y, ds_clim.pr, kind="linear")
    ds_clim = xr.open_dataset(os.path.join(ATMOS_DATA_PATH, "ps-ECMWF-clim.nc"))
    fsp = interp2d(ds_clim.X, ds_clim.Y, ds_clim.ps, kind="linear")

    wnsp = fwnsp(x_axis, y_axis_u)
    wnsp[wnsp < wnsp_min] = wnsp_min
    ds["wnspClim"] = (["Yu", "X"], wnsp)
    ds["tsClim"] = (["Yu", "X"], fts(x_axis, y_axis_u))
    ds["prClim"] = (["Yu", "X"], fpr(x_axis, y_axis_u))
    ds["spClim"] = (["Yu", "X"], fsp(x_axis, y_axis_u))

    # TRENDS
    ds_trend = xr.open_dataset(os.path.join(ATMOS_DATA_PATH, "ts-ECMWF-trend.nc"))
    fts_trend = interp2d(ds_trend.X, ds_trend.Y, ds_trend.ts, kind="linear")
    ds_trend = xr.open_dataset(os.path.join(ATMOS_DATA_PATH, "pr-ECMWF-trend.nc"))
    fpr_trend = interp2d(ds_trend.X, ds_trend.Y, ds_trend.pr, kind="linear")

    ts_trend = fts_trend(x_axis, y_axis_u)
    ds["tsTrend"] = (["Yu", "X"], ts_trend)

    pr_trend = fpr_trend(x_axis, y_axis_u)
    pr_trend[abs(y_axis_u) > 25] = 0
    pr_trend[pr_trend > 5e-5] = 5e-5
    ds["prTrend"] = (["Yu", "X"], pr_trend)
    ds["prTrend"] = smooth121(ds.prTrend, ["Yu", "X"], perdims=["X"])

    dsmask = xr.open_dataset(os.path.join(ATMOS_DATA_PATH, "mask-360x180.nc"))
    fmask = interp2d(dsmask.X, dsmask.Y, dsmask.mask, kind="linear")
    ds["mask"] = (["Yu", "X"], fmask(x_axis, y_axis_u))

    # tsClim = ds.tsClim.values
    sp_clim = ds.spClim.values
    wnsp_clim = ds.wnspClim.values
    wnsp_clim[wnsp_clim < wnsp_min] = wnsp_min
    mask = ds.mask.values
    wend = wnsp_clim
    wbeg = wnsp_clim

    tsend = (ds.tsClim + (1 - mask) * ds.tsTrend / 2).values
    ts_beg = (ds.tsClim - (1 - mask) * ds.tsTrend / 2).values
    prend = (ds.prClim + ds.prTrend / 2).values
    prbeg = (ds.prClim - ds.prTrend / 2).values
    q_th_end = K1 * (tsend - 30) / b_coeff
    q_th_beg = K1 * (ts_beg - 30) / b_coeff

    qaend = f_qa(tsend, sp_clim)
    # qaend = f_qa2(tsend)
    e_end = f_evap(mask, qaend, wnsp_clim)
    pr_end = e_end
    pr_end[pr_end < 0] = 0
    # pr_end[pr_end>pr_max] = pr_max

    qa_beg = f_qa(ts_beg, sp_clim)
    # qa_beg = f_qa2(ts_beg)
    e_beg = f_evap(mask, qa_beg, wnsp_clim)
    pr_beg = e_beg
    pr_beg[pr_beg < 0] = 0
    # pr_beg[pr_beg>pr_max] = pr_max

    q_th = q_th_end
    pr = pr_end
    e1 = e_end
    qa1 = qaend

    # Find total pr, u and v at end
    for _ in range(0, number_iterations):
        # Start main calculation
        q_c = (
            np.pi * latent_heat_vap * pr / (2 * cp_air * rho_00 * height_tropopause)
        )  # heating from precip
        # convective heating part, Qc
        q1 = b_coeff * (q_c + q_th)
        # Q1 is a modified heating since the part
        # involving θ is on the left-hand side
        (u1, v1, phi1) = s91_solver(q1)
        d_amc = xr.DataArray(f_mc(qa1, u1, v1), dims=["Yu", "X"])
        mc1 = smooth121(d_amc, ["Yu", "X"], perdims=["X"]).values
        if prcp_land:
            pr = (1 - mask) * (mc1 + e1) + mask * prend
        else:
            pr = (1 - mask) * (mc1 + e1)
        pr[pr < 0] = 0
        # pr[pr > pr_max] = pr_max

    mc_end = mc1
    uend = u1
    vend = v1
    phiend = phi1
    pr_end = pr

    q_th = q_th_beg
    pr = pr_beg
    e1 = e_beg
    qa1 = qa_beg

    # Find total pr, u and v at beginning
    for _ in range(0, number_iterations):
        # Start main calculation
        q_c = (
            np.pi * latent_heat_vap * pr / (2 * cp_air * rho_00 * height_tropopause)
        )  # heating from precip
        q1 = b_coeff * (q_c + q_th)
        (u1, v1, phi1) = s91_solver(q1)
        d_amc = xr.DataArray(f_mc(qa1, u1, v1), dims=["Yu", "X"])
        mc1 = smooth121(d_amc, ["Yu", "X"], perdims=["X"]).values
        if prcp_land:
            pr = (1 - mask) * (mc1 + e1) + mask * prbeg
        else:
            pr = (1 - mask) * (mc1 + e1)
        pr[pr < 0] = 0
        # pr[pr > pr_max] = pr_max

    mc_beg = mc1
    ubeg = u1
    vbeg = v1
    phibeg = phi1
    pr_beg = pr

    # save and plot the trends
    ds["utrend"] = (["Yu", "X"], uend - ubeg)
    ds["vtrend"] = (["Yv", "X"], vend - vbeg)
    ds["phitrend"] = (["Yu", "X"], phiend - phibeg)
    ds["tstrend"] = (["Yu", "X"], tsend - ts_beg)
    ds["PRtrend"] = (["Yu", "X"], pr_end - pr_beg)
    ds["Qthtrend"] = (["Yu", "X"], q_th_end - q_th_beg)
    ds["uend"] = (["Yu", "X"], uend)
    ds["vend"] = (["Yv", "X"], vend)
    ds["wend"] = (["Yu", "X"], wend)
    ds["phiend"] = (["Yu", "X"], phiend)
    ds["tsend"] = (["Yu", "X"], tsend)
    ds["PRend"] = (["Yu", "X"], pr_end)
    ds["Qthend"] = (["Yu", "X"], q_th_end)
    ds["Eend"] = (["Yu", "X"], e_end)
    ds["MCend"] = (["Yu", "X"], mc_end)
    ds["qaend"] = (["Yu", "X"], qaend)
    ds["ubeg"] = (["Yu", "X"], ubeg)
    ds["vbeg"] = (["Yv", "X"], vbeg)
    ds["wbeg"] = (["Yu", "X"], wbeg)
    ds["phibeg"] = (["Yu", "X"], phibeg)
    ds["tsbeg"] = (["Yu", "X"], ts_beg)
    ds["PRbeg"] = (["Yu", "X"], pr_beg)
    ds["Qthbeg"] = (["Yu", "X"], q_th_beg)
    ds["Ebeg"] = (["Yu", "X"], e_beg)
    ds["MCbeg"] = (["Yu", "X"], mc_beg)
    ds["qabeg"] = (["Yu", "X"], qa_beg)

    # There is 2 gridpoint noise in the phi field - so add a smooth in X:
    ds["phitrend"] = smooth121(ds.phitrend, ["X"], number_smooths=1, perdims=["X"])

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

    for diff_direc in [ATMOS_TMP_PATH, direc]:

        basedir = os.path.join(diff_direc, "S91")

        if diff_direc != "":
            if not os.path.isdir(diff_direc):
                os.makedirs(diff_direc)

        outfile = basedir + "-hq" + str(h_q) + "-prcp_land" + str(prcp_land) + ".nc"

        print(outfile)

        ds.to_netcdf(outfile, encoding=en_dict)


###--------------------------- Begin dQ ----------------------------


@timeit
def get_dclim(direc: str = "") -> any:
    """Opens the files, and applies functions.

    Returns:
        any: A list of outputs.
            dclim, u_b, alh, alw, blw, dtemp_se, rh, c_b, t_sb.
    """
    files = []

    for i, m in enumerate(mem):
        name = names[m]
        variable = var[i]
        file = os.path.join(ATMOS_DATA_PATH, variable + "-" + name + "-clim60.nc")
        print(name, variable, file)
        print(file)
        assert os.path.isfile(file)
        files += [file]

    dclim_loc = xr.open_mfdataset(files, decode_times=False)

    # set Q'_LW + Q'_LH = 0, solve for Ts' (assuming U'=0)
    # Q'_LW = (alw(temp_surface_bar,c_bar,f1_bar)* Tsprime
    #          + blw(temp_surface_bar,c_bar) * f1prime)
    # Q'_LH = alh(temp_surface_bar,u_bar) * Tsprime
    # Q'_LH is from formula 13 in paper
    # Q'_LW is from formula 14 in paper

    t_sb_loc = 1.0 * dclim_loc.ts
    tmp = 1.0 * dclim_loc.sfcWind.stack(z=("lon", "lat")).load()
    tmp[tmp < wnsp_min] = wnsp_min
    u_b_loc = tmp.unstack("z").T
    c_b_loc = dclim_loc.clt / 100.0
    rh_loc = dclim_loc.rh / 100.0
    f1p = -0.003

    alh0 = f_dqlh_dtemp(t_sb_loc, u_bar, rh_loc)
    alw0 = f_dqlw_dtemp(t_sb_loc, c_bar, f1_bar, rh_loc)
    blw0 = f_dqlw_df(t_sb_loc, c_bar)
    dtemp_se0 = -blw0 * f1p / (alh0 + alw0)
    dclim_loc["dTse0"] = dtemp_se0

    alh1 = f_dqlh_dtemp(t_sb_loc, u_b_loc, rh_loc)
    alw1 = f_dqlw_dtemp(t_sb_loc, c_bar, f1_bar, rh_loc)
    blw1 = f_dqlw_df(t_sb_loc, c_bar)
    dtemp_se1 = -blw1 * f1p / (alh1 + alw1)
    dclim_loc["dTse1"] = dtemp_se1

    alh2 = f_dqlh_dtemp(t_sb_loc, u_bar, rh_loc)
    alw2 = f_dqlw_dtemp(t_sb_loc, c_b_loc, f1_bar, rh_loc)
    blw2 = f_dqlw_df(t_sb_loc, c_b_loc)
    dtemp_se2 = -blw2 * f1p / (alh2 + alw2)
    dclim_loc["dTse2"] = dtemp_se2

    alh_loc = f_dqlh_dtemp(t_sb_loc, u_b_loc, rh_loc)
    alw_loc = f_dqlw_dtemp(t_sb_loc, c_b_loc, f1_bar, rh_loc)
    blw_loc = f_dqlw_df(t_sb_loc, c_b_loc)
    dtemp_se_loc = -blw_loc * f1p / (alh_loc + alw_loc)

    dclim_loc["dTse"] = dtemp_se_loc
    dclim_loc["ALH"] = alh_loc
    dclim_loc["ALW"] = alw_loc
    dclim_loc["BLW"] = blw_loc
    dclim_loc["QLW"] = alw_loc + blw_loc * f1p / dtemp_se_loc

    dclim_loc.to_netcdf(os.path.join(direc, "Q.nc"))

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


def make_figure(
    direc: str = "",
    cmap: Union[str] = "viridis",
    lat: str = "latitude",
    lon: str = "longitude",
) -> None:
    """Make figure.

    Args:
    Args:
        direc (str): direc.
        cmap (str, optional): matplotlib colormap. Defaults to "viridis".
        lat (str, optional): [description]. Defaults to "latitude".
        lon (str, optional): [description]. Defaults to "longitude".
    """
    log("Make figure.")

    dclim, u_b, _, _, _, _, _, _, _ = get_dclim(direc=direc)

    v_min_ad = 0.0
    v_max_ad = 0.6

    plt.figure(figsize=(8, 6))
    plt.subplot(321)
    dp = dclim.dTse0.plot.contourf(
        levels=11, cmap=cmap, vmin=v_min_ad, vmax=v_max_ad, add_colorbar=0
    )
    # ,vmin=-2,vmax=2,add_colorbar=0)
    plt.title(r"$T^{\,\prime}_s$  for $\bar U,\bar C$")
    plt.ylabel(lat)
    plt.xlabel(lon)
    # cbar = plt.colorbar(dp)
    _ = plt.colorbar(dp)
    plt.subplot(322)
    dp = dclim.dTse1.plot.contourf(
        levels=11, cmap=cmap, vmin=v_min_ad, vmax=v_max_ad, add_colorbar=0
    )
    # ,vmin=-2,vmax=2,add_colorbar=0)
    plt.title(r"$T^{\,\prime}_s$  for $\bar U(x,y), \bar C$")
    plt.ylabel(lat)
    plt.xlabel(lon)
    _ = plt.colorbar(dp)
    plt.subplot(323)
    dp = dclim.dTse2.plot.contourf(
        levels=11, cmap=cmap, vmin=v_min_ad, vmax=v_max_ad, add_colorbar=0
    )
    # ,vmin=-2,vmax=2,add_colorbar=0)
    plt.title(r"$T^{\,\prime}_s$  for $\bar U, \bar C(x,y)$")
    plt.ylabel(lat)
    plt.xlabel(lon)
    _ = plt.colorbar(dp)
    plt.subplot(324)
    dp = dclim.dTse.plot.contourf(
        levels=11, cmap=cmap, vmin=v_min_ad, vmax=v_max_ad, add_colorbar=0
    )
    # ,vmin=-2,vmax=2,add_colorbar=0)
    plt.title(r"$T^{\,\prime}_s$  for $\bar U(x,y), \bar C(x,y)$")
    plt.ylabel(lat)
    plt.xlabel(lon)
    _ = plt.colorbar(dp)
    plt.subplot(325)
    dp = (dclim.clt / 100).plot.contourf(
        levels=11, cmap=cmap, vmin=0.0, vmax=1.0, add_colorbar=0
    )
    # ,vmin=-2,vmax=2,add_colorbar=0)
    plt.title(r"$\bar C(x,y)$")
    plt.ylabel(lat)
    plt.xlabel(lon)
    _ = plt.colorbar(dp)
    plt.subplot(326)
    dp = u_b.plot.contourf(levels=11, cmap=cmap, vmin=4.0, vmax=8.0, add_colorbar=0)
    # ,vmin=-2,vmax=2,add_colorbar=0)
    plt.title(r"$\bar U(x,y)$")
    plt.ylabel(lat)
    plt.xlabel(lon)
    _ = plt.colorbar(dp)
    plt.tight_layout()
    plt.savefig(os.path.join(ATMOS_PATH, "Tsp4.eps"), format="eps", dpi=1000)
    # plt.show()


def output_dq(direc: str = "") -> None:
    """Outputs "dQ.nc".

    Args:
        direc ([type]): [description]
    """

    dclim, u_b, alh, alw, blw, dtemp_se, rh, c_b, t_sb = get_dclim(direc=direc)

    # Now, save the dq_df and dq_dt terms for using in TCOM:
    dq_dt = alh + alw
    dq_df = blw

    # Define the new Dataset
    dq = xr.Dataset(
        {
            "lon": ("lon", dclim.lon),
            "lat": ("lat", dclim.lat),
            "dq_dt": (["lat", "lon"], dq_dt),
            "dq_df": (["lat", "lon"], dq_df),
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
    dq.to_netcdf(os.path.join(direc, "dQ.nc"))


if __name__ == "__main__":
    # get_dclim()
    # python3 src/models/atmos.py
    output_trends(direc="")
    output_dq(direc="")
    make_figure(direc="")
