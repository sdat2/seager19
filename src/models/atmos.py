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
from typing import Tuple
import os
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from scipy.interpolate import interp2d
from scipy.fftpack import fft, ifft
from src.constants import ATMOS_TMP_PATH, ATMOS_DATA_PATH, ATMOS_PATH, PROJECT_PATH


# ------------- constants -----------------------
# begining TCAM
k_days = 10
efrac = 2.0  # multiply epsu by efrac to get epsv
hq = 1800  # Hq is a scale depth for moisture
prcp_land = 1  # use data precip trends over land
wnsp_min = 4
rho00 = 0.3  # this is actual rho_bar in the paper
# 0.3 kg m-3
prmax = 20.0 / 3600 / 24
r = 0.80  # relative humidity uniformly 0.8
number_iterations = 50
gravity = 9.8  #  m s-2
height_tropopause = 15000  # metres
theta_00 = 300
nbsq = 3.0e-4
radius_earth = 6.37e6  # metres
omega2 = 2 * (2 * np.pi / 86400)
latent_heat_vap = 2.5e6  # latent heat
cp_air = 1000  #  cp_air is the specific heat capacity of air.
b_coeff = gravity * np.pi / (nbsq * theta_00 * height_tropopause)
eps_days = 0.75
eps = 1.0 / (eps_days * 86400)  # 1/.75 d
epsu = eps  # 1/.75 d
epsv = efrac * eps  # efrac=1/2 in paper
K1 = b_coeff / (k_days * 86400)
epsp = (np.pi / height_tropopause) ** 2 / (nbsq * k_days * 86400)
beta = omega2 / radius_earth
rho_air = 1.225  # kg m-3 - also called rho_00
c_e = 0.00125
eps = 0.97  # problem here with second definintion.
stefan_boltzman_const = 5.67e-8
ps = 1000  # pressure at the surface?
es0 = 6.11
delta = 1.0
f2 = 0.05
# 'a' should decrease when deep convection happens above 28 degC
#  a = Ts-temp_0c;a[a>28] = 40;a[a<=28] = 80;a = 0.01*a
a = 0.6  # this isn't the option used in the paper.

# basic parameters
temp_0c = 273.15
f1_bar = 0.39
u_bar = 5.0
temp_surface_bar = temp_0c + 25
c_bar = 0.6

# grid characteristics

nx = 180
ny = 60
y_north_lim = 60
y_south_lim = -y_north_lim

# make grids
dx = 360 / nx
dy = (y_north_lim - y_south_lim) / ny
x_axis = np.linspace(0, 360 - dx, nx)
y_axis_v = np.linspace(y_south_lim + dy / 2, y_north_lim - dy / 2, ny)
y_axis_u = np.linspace(y_south_lim + dy, y_north_lim - dy, ny - 1)
y_axis_i = np.linspace(y_south_lim + 3 * dy / 2, y_north_lim - 3 * dy / 2, ny - 2)
x_spacing = x_axis[1] - x_axis[0]  # degrees
y_spacing = y_axis_v[1] - y_axis_v[0]  # degrees
dxm = x_spacing * radius_earth * np.pi / 180
dym = y_spacing * radius_earth * np.pi / 180
dym2 = dym * dym

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
const1 = rho_air * c_e * latent_heat_vap

mem = "EEEf"

# the different model names in a dict?
names = {
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

var = {0: "ts", 1: "clt", 2: "sfcWind", 3: "rh"}


# --------------- flux functions ----------------------


def f_cor(y: np.ndarray) -> np.ndarray:
    """Corriolis force coeff.

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
    return es0 * np.exp(
        17.67 * (temperature - temp_0c) / (temperature - temp_0c + 243.5)
    )


def f_qs(temperature: xr.DataArray) -> xr.DataArray:
    """flux qs.

    Args:
        temperature (xr.DataArray): temp.

    Returns:
        xr.DataArray: Flux qs.
    """
    return 0.622 * f_es(temperature) / ps


def f_dqs_dtemp(temperature: xr.DataArray) -> xr.DataArray:
    """flux dqs.

    Args:
        temperature (xr.DataArray): temp.

    Returns:
        xr.DataArray: flux qs.
    """
    return f_qs(temperature) * (17.67 * 243.5) / (temperature - temp_0c + 243.5) ** 2


def f_qlh(
    temperature: xr.DataArray, u_sp: xr.DataArray, rh_loc: xr.DataArray
) -> xr.DataArray:
    """flux qlh.

    Args:
        temperature (xr.DataArray): [description]
        u_sp (xr.DataArray): [description]
        rh_loc (xr.DataArray): [description]

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
        temperature (xr.DataArray): [description]
        u_sp (xr.DataArray): [description]
        rh_loc (xr.DataArray): [description]

    Returns:
        xr.DataArray: flux dqlh_dtemp.
    """
    # print("u_sp", type(u_sp))
    return const1 * u_sp * f_dqs_dtemp(temperature) * (1 - rh_loc)


# Find linearization of Q_LW (longwave)
const2 = eps * stefan_boltzman_const


def f_temp_a(temperature: xr.DataArray) -> xr.DataArray:
    """temperature anomaly.

    Args:
        temperature (xr.DataArray): temp.

    Returns:
        xr.DataArray: temperature anomaly.
    """
    return temperature - delta


def f_ebar(temperature: xr.DataArray, rh_loc: xr.DataArray) -> xr.DataArray:
    """flux e_bar.

    Args:
        temperature (xr.DataArray): [description]
        rh_loc (xr.DataArray): [description]

    Returns:
        xr.DataArray: [description]
    """
    qa = rh_loc * f_qs(temperature)
    return qa * ps / 0.622


def f_qlw1(
    temperature: xr.DataArray, const: xr.DataArray, f: float, rh_loc: xr.DataArray
) -> xr.DataArray:
    """[summary]

    Args:
        temperature (xr.DataArray): [description]
        const (xr.DataArray): [description]
        f (float): [description]
        rh_loc (xr.DataArray): [description]

    Returns:
        xr.DataArray: [description]
    """
    # print("rh_loc", type(rh_loc))
    temp_a = f_temp_a(temperature)
    return (
        const2
        * (1 - a * const ** 2)
        * temp_a ** 4
        * (f - f2 * np.sqrt(f_ebar(temperature, rh_loc)))
    )


def f_qlw2(temperature: xr.DataArray) -> xr.DataArray:
    """[summary]

    Args:
        temperature (xr.DataArray): [description]

    Returns:
        xr.DataArray: [description]
    """
    # print("temperature", type(temperature))
    # print("const", type(const))
    return (
        4
        * eps
        * stefan_boltzman_const
        * temperature ** 3
        * (temperature - f_temp_a(temperature))
    )


# def f_QLW(T, f, rh_loc):
#     return f_qlw1(T, f, rh_loc) + f_qlw2(T)
# This function seemed to call a function in an incorrect way.


def f_dqlw_df(temperature: xr.DataArray, const: xr.DataArray) -> xr.DataArray:
    """flux dqlw_df.

    Args:
        temperature (xr.DataArray): temp.
        const (xr.DataArray): constant.

    Returns:
        xr.DataArray: flux dqlw_df.
    """
    return const2 * (1 - a * const ** 2) * temperature ** 4


def f_dqlw_dtemp(
    temperature: xr.DataArray, const: xr.DataArray, f: float, rh_loc: xr.DataArray
) -> xr.DataArray:
    """flux dqlw_dtemp.

    Args:
        temperature (xr.DataArray): [description]
        const (xr.DataArray): [description]
        f (float): [description]
        rh_loc (xr.DataArray): [description]

    Returns:
        xr.DataArray: [description]
    """
    ebar = f_ebar(temperature, rh_loc)
    qs = f_qs(temperature)
    dqs_dtemp = f_dqs_dtemp(temperature)
    return const2 * (
        (1 - a * const ** 2)
        * temperature ** 3
        * (4 * f - f2 * np.sqrt(ebar) * (4 + temperature * dqs_dtemp / 2 / qs))
        + 12 * temperature ** 2 * delta
    )


def f_qa(ts: np.ndarray, sp: np.ndarray) -> np.ndarray:
    """f_qa.

    Args:
        ts (np.ndarray): sst in Kelvin
        sp (np.ndarray): surface pressure in mb

    Returns:
        np.ndarray: qs, surface specific humidity
    """
    efac = 0.622
    es = es0 * np.exp(17.67 * (ts - 273.15) / ((ts - 273.15) + 243.5))
    return efac * r * es / sp


def f_qa2(temp_surface: np.ndarray) -> np.ndarray:
    """flux qa2.

    Args:
        temp_surface (np.ndarray): sst in Kelvin.

    Returns:
        np.ndarray: qs, surface specific humidity.
    """
    return 0.001 * (temp_surface - 273.15 - 11.0)


def f_evap(mask: np.ndarray, qa: np.ndarray, wnsp: np.ndarray) -> np.ndarray:
    """evaporation flux.

    Args:
        mask (np.ndarray): [description]
        qa (np.ndarray): surface air humidity
        wnsp (np.ndarray): surface windspeed in m/s

    Returns:
        np.ndarray: Evap in kg/m^2/s
    """
    c_s_e = 0.0015 * (1 + mask / 2)
    # c_s_e = 0.0012
    return c_s_e * rho_air * (1 - r) * qa * wnsp / r


def f_mc(qa: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """moisture convergence flux.

    To calculate surface heat fluxes and atmospheric moisture
    convergence, relative humidity is assumed to be spatially
    uniform in our standard model.
    (N.B., v is on y_axis_v points, u,q are on y_axis_u points)

    Args:
        qa (np.ndarray): surface air humidity
        u (np.ndarray): low level winds in m/s
        v (np.ndarray): low level winds in m/s

    Returns:
        np.ndarray: Moisture Convergence in kg/m^2/s
    """
    qu = qa * u
    qux = ifft(1.0j * kk_wavenumber * fft(qu) / radius_earth).real
    aq = (qa[1 : ny - 1, :] + qa[0 : ny - 2, :]) / 2.0
    qv = aq * v[1 : ny - 1, :]
    z = np.zeros((1, nx))
    qv = np.concatenate((z, qv, z), axis=0)
    # qvy = qv.diff('Yu')/dym
    qvy = (qv[1:ny, :] - qv[0 : ny - 1, :]) / dym
    return -hq * (qux + qvy) * rho_air


# ---------------- equation solvers ---------------------


def tdma_solver(
    ny_loc: int, a_loc: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """tdma solver.

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


def s91_solver(q1: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """S91 folder from TCAM.py.


    The atmosphere equations are solved by Fourier transforming in longitude,
    forming an equation for v for each zonal wavenumber that is finite
    differenced, and the resulting tri-diagonal system is solved by matrix
    inversion, transforming back into longitude.

    Usef fft, ifft.

               g pi N ^ 2
        q1 = --------------(k theta_s Q_c)
              theta_00 z_t

    Args:
        q1 (np.ndarray): modified heating that drives winds.

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: (u, v, phi)
    """
    # print("q1", type(q1))
    q1_time = fft(q1)
    f_q = fcu[:, np.newaxis] * q1_time
    a_f_q = (f_q[1 : ny - 1, :] + f_q[0 : ny - 2, :]) / 2.0
    km = kk_wavenumber / radius_earth
    d_q = (q1_time[1 : ny - 1, :] - q1_time[0 : ny - 2, :]) / dym
    rk = 1.0j * km * beta - epsu * epsv * epsp - epsv * km ** 2
    fcp = fcu[1 : ny - 1] ** 2 / 4.0
    fcm = fcu[0 : ny - 2] ** 2 / 4.0

    ak = epsu / dym2 - epsp * fcm[:, np.newaxis]
    ck = epsu / dym2 - epsp * fcp[:, np.newaxis]
    bk = (
        -2 * epsu / dym2
        - epsp * (fcm[:, np.newaxis] + fcp[:, np.newaxis])
        + rk[np.newaxis, :]
    )
    dk = -epsu * d_q + 1.0j * km[np.newaxis, :] * a_f_q

    vtk = tdma_solver(ny - 2, ak, bk, ck, dk)

    z = np.zeros((1, nx))
    vt = np.concatenate((z, vtk, z), axis=0)
    av = (vt[1:ny, :] + vt[0 : ny - 1, :]) / 2.0
    fav = fcu[:, np.newaxis] * av
    dv = (vt[1:ny, :] - vt[0 : ny - 1, :]) / dym
    coeff = epsu * epsp + km * km
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


def output_trends() -> None:
    """output trends ds"""
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
    ds["epsv"] = eps_days / efrac
    ds.epsv.attrs = [("units", "day")]
    ds["hq"] = hq
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
    # pr_end[pr_end>prmax] = prmax

    qa_beg = f_qa(ts_beg, sp_clim)
    # qa_beg = f_qa2(ts_beg)
    e_beg = f_evap(mask, qa_beg, wnsp_clim)
    pr_beg = e_beg
    pr_beg[pr_beg < 0] = 0
    # pr_beg[pr_beg>prmax] = prmax

    q_th = q_th_end
    pr = pr_end
    e1 = e_end
    qa1 = qaend

    # Find total pr, u and v at end
    for _ in range(0, number_iterations):
        # Start main calculation
        q_c = (
            np.pi * latent_heat_vap * pr / (2 * cp_air * rho00 * height_tropopause)
        )  # heating from precip
        q1 = b_coeff * (q_c + q_th)
        (u1, v1, phi1) = s91_solver(q1)
        d_amc = xr.DataArray(f_mc(qa1, u1, v1), dims=["Yu", "X"])
        mc1 = smooth121(d_amc, ["Yu", "X"], perdims=["X"]).values
        if prcp_land:
            pr = (1 - mask) * (mc1 + e1) + mask * prend
        else:
            pr = (1 - mask) * (mc1 + e1)
        pr[pr < 0] = 0
        # pr[pr > prmax] = prmax

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
            np.pi * latent_heat_vap * pr / (2 * cp_air * rho00 * height_tropopause)
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
        # pr[pr > prmax] = prmax

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

    basedir = os.path.join(ATMOS_TMP_PATH, "S91")

    if not os.path.isdir(ATMOS_TMP_PATH):
        os.makedirs(ATMOS_TMP_PATH)

    outfile = basedir + "-hq" + str(hq) + "-prcp_land" + str(prcp_land) + ".nc"
    print(outfile)

    en_dict = {
        "K": {"dtype": "f4"},
        "epsu": {"dtype": "f4"},
        "epsv": {"dtype": "f4"},
        "hq": {"dtype": "f4"},
    }

    ds.to_netcdf(outfile, encoding=en_dict)


###--------------------------- Begin dQ ----------------------------


def get_dclim() -> any:
    """opens the files, and applies functions

    Returns:
        any: great.
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

    dclim_loc.to_netcdf(os.path.join(PROJECT_PATH, "Q.nc"))

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


dclim, u_b, alh, alw, blw, dtemp_se, rh, c_b, t_sb = get_dclim()


def make_figure(
    cmap: str = "viridis", lat: str = "latitude", lon: str = "longitude"
) -> None:
    """Make figure.

    Args:
        cmap (str, optional): matplotlib colormap. Defaults to "viridis".

    """
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


def output_dq() -> None:
    """outputs "dQ.nc"."""
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
    dq.to_netcdf(os.path.join(PROJECT_PATH, "dQ.nc"))


if __name__ == "__main__":
    # get_dclim()
    output_trends()
    output_dq()
    make_figure()
