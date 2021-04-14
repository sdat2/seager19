"""src.models.atmos.

pytest src/test/test_atmos.py
python3 src/models/atmos.py
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
eps_days = 0.75
k_days = 10
efrac = 2.0  # multiply epsu by efrac to get epsv
hq = 1800
prcp_land = 1  # use data precip trends over land
wnspmin = 4
rho00 = 0.3
prmax = 20.0 / 3600 / 24
r = 0.80
nx = 180
ny = 60
YN = 60
YS = -YN
number_iterations = 50
gravity = 9.8
height_tropopause = 15000
theta_00 = 300
NBSQ = 3.0e-4
radius_earth = 6.37e6  # m
omega2 = 2 * (2 * np.pi / 86400)
L = 2.5e6  # latent heat
cp_air = 1000
B = gravity * np.pi / (NBSQ * theta_00 * height_tropopause)
eps = 1.0 / (eps_days * 86400)
epsu = eps
epsv = efrac * eps
K1 = B / (k_days * 86400)
epsp = (np.pi / height_tropopause) ** 2 / (NBSQ * k_days * 86400)
beta = omega2 / radius_earth

# make grids
dx = 360 / nx
dy = (YN - YS) / ny
y_axis_v = np.linspace(YS + dy / 2, YN - dy / 2, ny)
x_axis = np.linspace(0, 360 - dx, nx)
y_axis_u = np.linspace(YS + dy, YN - dy, ny - 1)
Yi = np.linspace(YS + 3 * dy / 2, YN - 3 * dy / 2, ny - 2)
x_spacing = x_axis[1] - x_axis[0]  # degrees
y_spacing = y_axis_v[1] - y_axis_v[0]  # degrees
dxm = x_spacing * radius_earth * np.pi / 180
dym = y_spacing * radius_earth * np.pi / 180
dym2 = dym * dym

# need to have the correct ordering of the wave numbers for fft
N = nx
if N % 2 == 0:
    Kk = np.asarray(
        list(range(0, N // 2)) + [0] + list(range(-N // 2 + 1, 0)), np.float64
    )
else:
    Kk = np.asarray(
        list(range(0, (N - 1) // 2)) + [0] + list(range(-(N - 1) // 2, 0)), np.float64
    )

rho_air = 1.225
c_e = 0.00125
eps = 0.97
stefan_boltzman_const = 5.67e-8
ps = 1000  # pressure at the surface?
es0 = 6.11
delta = 1.0
f2 = 0.05
# 'a' should decrease when deep convection happens above 28 degC
# a = Ts-temp_0c;a[a>28] = 40;a[a<=28] = 80;a = 0.01*a
a = 0.6

# basic parameters
temp_0c = 273.15
f1bar = 0.39
u_bar = 5.0
temp_surface_bar = temp_0c + 25
c_bar = 0.6
wnspmin = 4.0

# Find linearization of Q_LH (latent heating)
const1 = rho_air * c_e * L

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
def fcor(y: np.ndarray) -> np.ndarray:
    """Corriolis force coeff.

    Args:
        y (np.ndarray): latitude

    Returns:
        np.ndarray: Corriolis force coeff.
    """
    return omega2 * y * np.pi / 180


fcu = fcor(y_axis_u)



def f_es(temperature: xr.DataArray) -> xr.DataArray:
    return es0 * np.exp(
        17.67 * (temperature - temp_0c) / (temperature - temp_0c + 243.5)
    )


def f_qs(temperature: xr.DataArray) -> xr.DataArray:
    return 0.622 * f_es(temperature) / ps


def f_dqsdT(temperature: xr.DataArray) -> xr.DataArray:
    return f_qs(temperature) * (17.67 * 243.5) / (temperature - temp_0c + 243.5) ** 2


def f_QLH(temperature: xr.DataArray, U: xr.DataArray, rh: xr.DataArray) -> xr.DataArray:
    # print("U", type(U))
    return const1 * U * f_qs(temperature) * (1 - rh)


def f_dQLHdT(
    temperature: xr.DataArray, U: xr.DataArray, rh: xr.DataArray
) -> xr.DataArray:
    # print("U", type(U))
    return const1 * U * f_dqsdT(temperature) * (1 - rh)


# Find linearization of Q_LW (longwave)
const2 = eps * stefan_boltzman_const


def f_Ta(temperature: xr.DataArray) -> xr.DataArray:
    return temperature - delta


def f_ebar(temperature: xr.DataArray, rh: xr.DataArray) -> xr.DataArray:
    qa = rh * f_qs(temperature)
    return qa * ps / 0.622


def f_QLW1(
    temperature: xr.DataArray, C: xr.DataArray, f: float, rh: xr.DataArray
) -> xr.DataArray:
    # print("rh", type(rh))
    Ta = f_Ta(temperature)
    return (
        const2
        * (1 - a * C ** 2)
        * Ta ** 4
        * (f - f2 * np.sqrt(f_ebar(temperature, rh)))
    )


def f_QLW2(temperature: xr.DataArray) -> xr.DataArray:
    # print("temperature", type(temperature))
    # print("C", type(C))
    return (
        4
        * eps
        * stefan_boltzman_const
        * temperature ** 3
        * (temperature - f_Ta(temperature))
    )


# def f_QLW(T, f, rh):
#     return f_QLW1(T, f, rh) + f_QLW2(T)


def f_dQLWdf(temperature: xr.DataArray, C: xr.DataArray) -> xr.DataArray:
    return const2 * (1 - a * C ** 2) * temperature ** 4


def f_dQLWdT(
    temperature: xr.DataArray, C: xr.DataArray, f: float, rh: xr.DataArray
) -> xr.DataArray:
    """[summary]

    Args:
        temperature (xr.DataArray): [description]
        C (xr.DataArray): [description]
        f (float): [description]
        rh (xr.DataArray): [description]

    Returns:
        xr.DataArray: [description]
    """
    ebar = f_ebar(temperature, rh)
    qs = f_qs(temperature)
    dqsdT = f_dqsdT(temperature)
    return const2 * (
        (1 - a * C ** 2)
        * temperature ** 3
        * (4 * f - f2 * np.sqrt(ebar) * (4 + temperature * dqsdT / 2 / qs))
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
    # ts: sst in Kelvin
    # return qs: surface specific humidity
    return 0.001 * (temp_surface - 273.15 - 11.0)


def f_E(mask: np.ndarray, qa: np.ndarray, wnsp: np.ndarray) -> np.ndarray:
    # qa: surface air humidity
    # wnsp: surface windspeed in m/s
    # return Evap in kg/m^2/s
    CsE = 0.0015 * (1 + mask / 2)
    # CsE = 0.0012
    return CsE * rho_air * (1 - r) * qa * wnsp / r


def f_MC(qa: np.ndarray, u: np.ndarray, v: np.ndarray) -> np.ndarray:
    # qa: surface air humidity
    # u,v: low level winds in m/s
    # (N.B., v is on y_axis_v points, u,q are on y_axis_u points)
    # return Moisture Convergence in kg/m^2/s
    qu = qa * u
    qux = ifft(1.0j * Kk * fft(qu) / radius_earth).real
    Aq = (qa[1 : ny - 1, :] + qa[0 : ny - 2, :]) / 2.0
    qv = Aq * v[1 : ny - 1, :]
    z = np.zeros((1, nx))
    qv = np.concatenate((z, qv, z), axis=0)
    # qvy = qv.diff('Yu')/dym
    qvy = (qv[1:ny, :] - qv[0 : ny - 1, :]) / dym
    return -hq * (qux + qvy) * rho_air


# ---------------- equation solvers ---------------------

def tdma_solver(
    nx: int, ny: int, a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """tdma solver

    Args:
        nx (int): number of equations
        ny (int): number of equations
        a (np.ndarray): [description]
        b (np.ndarray): [description]
        c (np.ndarray): [description]
        d (np.ndarray): [description]

    Returns:
        np.ndarray: xc
    """

    nf = ny  # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))  # copy arrays

    for it in range(1, nf):
        mc = ac[it, :] / bc[it - 1, :]
        bc[it, :] = bc[it, :] - mc * cc[it - 1, :]
        dc[it, :] = dc[it, :] - mc * dc[it - 1, :]

    xc = bc
    xc[-1, :] = dc[-1, :] / bc[-1, :]

    for il in range(nf - 2, -1, -1):
        xc[il, :] = (dc[il, :] - cc[il, :] * xc[il + 1, :]) / bc[il, :]

    return xc


def S91_solver(Q1: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """S91 folder from TCAM.py.

    Usef fft, ifft.

    Args:
        Q1 (np.ndarray): [description]

    Returns:
        Tuple[np.ndarray, np.ndarray, np.ndarray]: (u, v, phi)
    """
    # print("Q1", type(Q1))
    Q1t = fft(Q1)
    fQ = fcu[:, np.newaxis] * Q1t
    AfQ = (fQ[1 : ny - 1, :] + fQ[0 : ny - 2, :]) / 2.0
    km = Kk / radius_earth
    DQ = (Q1t[1 : ny - 1, :] - Q1t[0 : ny - 2, :]) / dym
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
    dk = -epsu * DQ + 1.0j * km[np.newaxis, :] * AfQ

    vtk = tdma_solver(nx, ny - 2, ak, bk, ck, dk)

    z = np.zeros((1, nx))
    vt = np.concatenate((z, vtk, z), axis=0)
    av = (vt[1:ny, :] + vt[0 : ny - 1, :]) / 2.0
    fav = fcu[:, np.newaxis] * av
    dv = (vt[1:ny, :] - vt[0 : ny - 1, :]) / dym
    coeff = epsu * epsp + km * km
    ut = (epsp * fav + 1.0j * (Q1t + dv) * km[np.newaxis, :]) / coeff[np.newaxis, :]
    phit = -(Q1t + 1.0j * ut * km[np.newaxis, :] + dv) / epsp
    v = ifft(vt).real
    u = ifft(ut).real
    phi = ifft(phit).real
    return (u, v, phi)


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
    wnsp[wnsp < wnspmin] = wnspmin
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
    spClim = ds.spClim.values
    wnspClim = ds.wnspClim.values
    wnspClim[wnspClim < wnspmin] = wnspmin
    mask = ds.mask.values
    wend = wnspClim
    wbeg = wnspClim

    tsend = (ds.tsClim + (1 - mask) * ds.tsTrend / 2).values
    tsbeg = (ds.tsClim - (1 - mask) * ds.tsTrend / 2).values
    prend = (ds.prClim + ds.prTrend / 2).values
    prbeg = (ds.prClim - ds.prTrend / 2).values
    q_th_end = K1 * (tsend - 30) / B
    q_th_beg = K1 * (tsbeg - 30) / B

    qaend = f_qa(tsend, spClim)
    # qaend = f_qa2(tsend)
    e_end = f_E(mask, qaend, wnspClim)
    pr_end = e_end
    pr_end[pr_end < 0] = 0
    # pr_end[pr_end>prmax] = prmax

    qabeg = f_qa(tsbeg, spClim)
    # qabeg = f_qa2(tsbeg)
    Ebeg = f_E(mask, qabeg, wnspClim)
    pr_beg = Ebeg
    pr_beg[pr_beg < 0] = 0
    # pr_beg[pr_beg>prmax] = prmax

    Qth = q_th_end
    PR = pr_end
    E1 = e_end
    qa1 = qaend

    # Find total PR, u and v at end
    for _ in range(0, number_iterations):
        # Start main calculation
        Qc = (
            np.pi * L * PR / (2 * cp_air * rho00 * height_tropopause)
        )  # heating from precip
        Q1 = B * (Qc + Qth)
        (u1, v1, phi1) = S91_solver(Q1)
        daMC = xr.DataArray(f_MC(qa1, u1, v1), dims=["Yu", "X"])
        MC1 = smooth121(daMC, ["Yu", "X"], perdims=["X"]).values
        if prcp_land:
            PR = (1 - mask) * (MC1 + E1) + mask * prend
        else:
            PR = (1 - mask) * (MC1 + E1)
        PR[PR < 0] = 0
        # PR[PR > prmax] = prmax

    MCend = MC1
    uend = u1
    vend = v1
    phiend = phi1
    pr_end = PR

    Qth = q_th_beg
    PR = pr_beg
    E1 = Ebeg
    qa1 = qabeg

    # Find total PR, u and v at beginning
    for _ in range(0, number_iterations):
        # Start main calculation
        Qc = (
            np.pi * L * PR / (2 * cp_air * rho00 * height_tropopause)
        )  # heating from precip
        Q1 = B * (Qc + Qth)
        (u1, v1, phi1) = S91_solver(Q1)
        daMC = xr.DataArray(f_MC(qa1, u1, v1), dims=["Yu", "X"])
        MC1 = smooth121(daMC, ["Yu", "X"], perdims=["X"]).values
        if prcp_land:
            PR = (1 - mask) * (MC1 + E1) + mask * prbeg
        else:
            PR = (1 - mask) * (MC1 + E1)
        PR[PR < 0] = 0
        # PR[PR > prmax] = prmax

    mc_beg = MC1
    ubeg = u1
    vbeg = v1
    phibeg = phi1
    pr_beg = PR

    # save and plot the trends
    ds["utrend"] = (["Yu", "X"], uend - ubeg)
    ds["vtrend"] = (["Yv", "X"], vend - vbeg)
    ds["phitrend"] = (["Yu", "X"], phiend - phibeg)
    ds["tstrend"] = (["Yu", "X"], tsend - tsbeg)
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
    ds["MCend"] = (["Yu", "X"], MCend)
    ds["qaend"] = (["Yu", "X"], qaend)
    ds["ubeg"] = (["Yu", "X"], ubeg)
    ds["vbeg"] = (["Yv", "X"], vbeg)
    ds["wbeg"] = (["Yu", "X"], wbeg)
    ds["phibeg"] = (["Yu", "X"], phibeg)
    ds["tsbeg"] = (["Yu", "X"], tsbeg)
    ds["PRbeg"] = (["Yu", "X"], pr_beg)
    ds["Qthbeg"] = (["Yu", "X"], q_th_beg)
    ds["Ebeg"] = (["Yu", "X"], Ebeg)
    ds["MCbeg"] = (["Yu", "X"], mc_beg)
    ds["qabeg"] = (["Yu", "X"], qabeg)

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

    dclim = xr.open_mfdataset(files, decode_times=False)

    # set Q'_LW + Q'_LH = 0, solve for Ts' (assuming U'=0)
    # Q'_LW = (ALW(temp_surface_bar,c_bar,f1bar)* Tsprime
    #          + BLW(temp_surface_bar,c_bar) * f1prime)
    # Q'_LH = ALH(temp_surface_bar,u_bar) * Tsprime

    t_sb = 1.0 * dclim.ts
    tmp = 1.0 * dclim.sfcWind.stack(z=("lon", "lat")).load()
    tmp[tmp < wnspmin] = wnspmin
    u_b = tmp.unstack("z").T
    c_b = dclim.clt / 100.0
    rh = dclim.rh / 100.0
    f1p = -0.003

    ALH0 = f_dQLHdT(t_sb, u_bar, rh)
    ALW0 = f_dQLWdT(t_sb, c_bar, f1bar, rh)
    BLW0 = f_dQLWdf(t_sb, c_bar)
    dTse0 = -BLW0 * f1p / (ALH0 + ALW0)
    dclim["dTse0"] = dTse0

    ALH1 = f_dQLHdT(t_sb, u_b, rh)
    ALW1 = f_dQLWdT(t_sb, c_bar, f1bar, rh)
    BLW1 = f_dQLWdf(t_sb, c_bar)
    dTse1 = -BLW1 * f1p / (ALH1 + ALW1)
    dclim["dTse1"] = dTse1

    ALH2 = f_dQLHdT(t_sb, u_bar, rh)
    ALW2 = f_dQLWdT(t_sb, c_b, f1bar, rh)
    BLW2 = f_dQLWdf(t_sb, c_b)
    dTse2 = -BLW2 * f1p / (ALH2 + ALW2)
    dclim["dTse2"] = dTse2

    ALH = f_dQLHdT(t_sb, u_b, rh)
    ALW = f_dQLWdT(t_sb, c_b, f1bar, rh)
    BLW = f_dQLWdf(t_sb, c_b)
    dTse = -BLW * f1p / (ALH + ALW)

    dclim["dTse"] = dTse
    dclim["ALH"] = ALH
    dclim["ALW"] = ALW
    dclim["BLW"] = BLW
    dclim["QLW"] = ALW + BLW * f1p / dTse
    dclim.to_netcdf(os.path.join(PROJECT_PATH, "Q.nc"))
    return dclim, u_b, ALH, ALW, BLW, dTse, rh, c_b, t_sb


dclim, u_b, ALH, ALW, BLW, dTse, rh, c_b, t_sb = get_dclim()


def make_figure(
    cmap: str = "viridis", lat: str = "latitude", lon: str = "longitude"
) -> None:
    """Make figure.

    Args:
        cmap (str, optional): matplotlib colormap. Defaults to "viridis".

    """
    plt.figure(figsize=(8, 6))
    plt.subplot(321)
    dp = dclim.dTse0.plot.contourf(
        levels=11, cmap=cmap, vmin=0.0, vmax=0.6, add_colorbar=0
    )
    # ,vmin=-2,vmax=2,add_colorbar=0)
    plt.title(r"$T^{\,\prime}_s$  for $\bar U,\bar C$")
    plt.ylabel(lat)
    plt.xlabel(lon)
    # cbar = plt.colorbar(dp)
    _ = plt.colorbar(dp)
    plt.subplot(322)
    dp = dclim.dTse1.plot.contourf(
        levels=11, cmap=cmap, vmin=0.0, vmax=0.6, add_colorbar=0
    )
    # ,vmin=-2,vmax=2,add_colorbar=0)
    plt.title(r"$T^{\,\prime}_s$  for $\bar U(x,y), \bar C$")
    plt.ylabel(lat)
    plt.xlabel(lon)
    _ = plt.colorbar(dp)
    plt.subplot(323)
    dp = dclim.dTse2.plot.contourf(
        levels=11, cmap=cmap, vmin=0.0, vmax=0.6, add_colorbar=0
    )
    # ,vmin=-2,vmax=2,add_colorbar=0)
    plt.title(r"$T^{\,\prime}_s$  for $\bar U, \bar C(x,y)$")
    plt.ylabel(lat)
    plt.xlabel(lon)
    _ = plt.colorbar(dp)
    plt.subplot(324)
    dp = dclim.dTse.plot.contourf(
        levels=11, cmap=cmap, vmin=0.0, vmax=0.6, add_colorbar=0
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
    dq_dt = ALH + ALW
    dq_df = BLW

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
    dq["ALH"] = ALH
    dq["ALW"] = ALW
    dq["BLW"] = BLW
    dq["dTse"] = dTse
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
