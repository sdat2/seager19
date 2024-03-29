"""Coupling between ocean and atmospheric models.

Example:

    Import statement usage::

        from src.models.coupling import Coupling

"""
from typing import Tuple, Union
from scipy.interpolate import interp2d
from scipy.constants import zero_Celsius
import xarray as xr
from typeguard import typechecked
from omegaconf import DictConfig
import wandb
from src.models.model_setup import ModelSetup
from src.models.atmos import Atmos
from src.models.ocean import Ocean
from src.visualisation.nino import get_nino_trend
from src.metrics import get_other_trends
from src.xr_utils import can_coords, open_dataset, cut_and_taper, get_trend
from src.visualisation.ani import animate_coupling
from src.visualisation.quiver import prcp_quiver_plot
from src.visualisation.trends import up_therm_qnet
from src.visualisation.comp_v_seager19 import (
    comp_oc_sst,
    comp_atm_prwnd,
    comp_oc_htherm,
)


# pylint: disable=no-value-for-parameter
class Coupling:
    """
    Coupled model part.

    Change over 1958–2017
    is computed by a linear trend. The atmosphere model is solved forced by a Ts
    comprised of the climatological mean for 1958–2017 plus and minus half of the
    SST trend and the difference of the two simulations taken to derive the change.
    For the coupled model, the ocean model is first forced with the change in CO2
    and climatological wind stress over 1958–2017. The resulting SST trend, plus
    the imposed heating change over land, are used to force the atmosphere model.
    The ocean model is forced again with both the changed wind stress and the CO2
    increase to derive a new SST change over 1958–2017 that is then used to force
    the atmosphere model. This iterative coupling is repeated until equilibrium is
    reached, which takes just a few times. There is a unique solution for any given
    value of CO2. The model wind stress change is computed as ρcWu*cD, where cD
    is a drag coefficient and u is the vector surface wind change computed by the
    atmosphere model, which is added to the ECMWF climatological stresses. Since
    the atmosphere model dynamics are only applicable in the tropics, the computed
    wind stress anomaly is only applied to the ocean model between 20S and 20N,
    and is linearly tapered to zero at 25S and 25N.
    """

    # pylint: disable=no-value-for-parameter
    def __init__(self, cfg: DictConfig, setup: ModelSetup) -> None:
        """Initialise model in standard way.

        Args:
            cfg (DictConfig): The config file for this model
                run containing run parameters.
            setup (ModelSetup): The setup object for
                this run containing parameters.

        """
        self.coup = cfg.coup
        self.cfg = cfg
        self.setup = setup
        self.ocean = Ocean(cfg, setup)
        self.atmos = Atmos(cfg, setup)

    @typechecked
    def f_stress(
        self,
        wind_speed_mean: Union[xr.DataArray, float],
        u_wind: xr.DataArray,
        v_wind: xr.DataArray,
    ) -> Tuple[xr.DataArray, xr.DataArray]:
        """Wind stress flux.

        .. math::
            :nowrap:

            \\begin{equation}
                \\vec{\\tau}= \\rho c_{D} W \\vec{u}
            \\end{equation}

        Args:
            wind_speed_mean (Union[xr.DataArray, float]): W, the climatological annual
                mean wind speed, which is taken from ECMWF
                reanalysis for our standard model and from the CMIP5
                multimodel mean when examining causes of bias in the
                CMIP5 model.
            u_wind (xr.DataArray): u_u, the zonal wind velocity.
            v_wind (xr.DataArray): u_v, the meridional wind velocity.

        Returns:
            Tuple[xr.DataArray, xr.DataArray]: tau_x (zonal wind stress),
                                               tau_y (meridional wind stress).

        """
        stress_coeff = self.coup.rho_air * self.coup.c_d * wind_speed_mean
        return stress_coeff * u_wind, stress_coeff * v_wind

    def get_tau_anom(
        self, wind: xr.DataArray, u_vel: xr.DataArray, v_vel: xr.DataArray
    ) -> Tuple[xr.DataArray, xr.DataArray]:
        """
        Return the tau anomaly with a clipping.

        Args:
            wind (xr.DataArray): wind speed field.
            u_vel (xr.DataArray): u wind velocity. (X, Yu).
            v_vel (xr.DataArray): v wind velocity. (X, Yv).

        Returns:
            Tuple[xr.DataArray, xr.DataArray]: tau_u, tau_v
        """
        sfcw50 = wind.sel(Y=slice(-50, 50))
        ds = xr.Dataset(
            {
                "X": ("X", sfcw50.X.values),
                "Y": ("Y", sfcw50.Y.values),
            }
        )
        fuend = interp2d(u_vel.X, u_vel.Yu, u_vel, kind="linear")
        ds["u_vel"] = (["Y", "X"], fuend(sfcw50.X.values, sfcw50.Y.values))
        fvend = interp2d(v_vel.X, v_vel.Yv, v_vel, kind="linear")
        ds["v_vel"] = (["Y", "X"], fvend(sfcw50.X.values, sfcw50.Y.values))
        t_u, t_v = self.f_stress(
            sfcw50,
            ds.u_vel,
            ds.v_vel,
        )
        return cut_and_taper(t_u).rename("tau_u"), cut_and_taper(t_v).rename("tau_v")

    def tau_anom_ds(self) -> xr.Dataset:
        """
        Wind stress anomaly.
        # TODO: Could simplify this function, and the following function.

        Returns:
            xr.Dataset: dataset with different different tau fields.

        """
        # var_list = [("ubeg", "vbeg"), ("utrend", "vtrend"), ("uend", "vend")]
        # new_var_list = [
        #    ("t_ubeg", "t_vbeg"),
        #    ("t_utrend", "t_vtrend"),
        #    ("t_uend", "t_vend"),
        # ]
        # ds = xr.open_dataset(self.setup.tcam_output())
        sfcwind = xr.open_dataset(self.setup.clim_file("sfcWind")).sfcWind
        ubeg = xr.open_dataset(self.setup.tcam_output()).ubeg
        vbeg = xr.open_dataset(self.setup.tcam_output()).vbeg
        utrend = xr.open_dataset(self.setup.tcam_output()).utrend
        vtrend = xr.open_dataset(self.setup.tcam_output()).vtrend
        uend = xr.open_dataset(self.setup.tcam_output()).uend
        vend = xr.open_dataset(self.setup.tcam_output()).vend
        t_beg_u, t_beg_v = self.get_tau_anom(sfcwind, ubeg, vbeg)
        t_end_u, t_end_v = self.get_tau_anom(sfcwind, uend, vend)
        t_trend_u, t_trend_v = self.get_tau_anom(sfcwind, utrend, vtrend)
        t_beg_u = t_beg_u.rename("t_beg_u")
        t_beg_v = t_beg_v.rename("t_beg_v")
        t_end_u = t_end_u.rename("t_end_u")
        t_end_v = t_end_u.rename("t_end_v")
        t_trend_u = t_trend_u.rename("t_trend_u")
        t_trend_v = t_trend_v.rename("t_trend_v")
        return xr.merge([t_beg_u, t_beg_v, t_end_u, t_end_v, t_trend_u, t_trend_v])

    def replace_stress(self, it: int) -> None:
        """Replace the stress files.

        Currently just resaves the clim files with a diff name.

        TODO: systematically explore other options:
            - linear change between taubeg and tauend
            - take half trend away or not.
            - turn off land precipitation.
            - make deep convection universal.
            - change in thermocline opposite of that expected - is this  a sign error?

        Args:
            it: the iteration in the coupling scheme.

        """
        # tau
        ds = self.tau_anom_ds()
        taux = xr.open_dataset(self.setup.tau_x(0), decode_times=False)
        taux_trend = ds.t_trend_u
        taux_beg = ds.t_beg_u
        taux_new = taux.copy()
        tauy = xr.open_dataset(self.setup.tau_y(0), decode_times=False)
        tauy_trend = ds.t_trend_v
        tauy_beg = ds.t_beg_v
        tauy_new = tauy.copy()
        # ah ok, this is definitely wrong
        time_length = len(tauy.coords["T"].values)
        # TODO: There is definitely a more elegant way of doing this.
        for i in range(time_length):
            if self.cfg.coup.add_stress:
                if self.cfg.coup.stress_trend:
                    taux_new["taux"][i, 0, 40:141, :] = (
                        taux.taux[i, 0, 40:141, :]
                        + (i / time_length - 1 / 2) * taux_trend[:, :]
                    )
                    tauy_new["tauy"][i, 0, 40:141, :] = (
                        tauy.tauy[i, 0, 40:141, :]
                        + (i / time_length - 1 / 2) * tauy_trend[:, :]
                    )
                else:
                    taux_new["taux"][i, 0, 40:141, :] = (
                        taux.taux[i, 0, 40:141, :]
                        + (i / time_length) * taux_trend[:, :]
                        + taux_beg[:, :]
                    )
                    tauy_new["tauy"][i, 0, 40:141, :] = (
                        tauy.tauy[i, 0, 40:141, :]
                        + (i / time_length) * tauy_trend[:, :]
                        + tauy_beg[:, :]
                    )

            else:
                if self.cfg.coup.stress_trend:
                    taux_new["taux"][i, 0, 40:141, :] = (
                        +(i / time_length - 1 / 2) * taux_trend[:, :]
                    )
                    tauy_new["tauy"][i, 0, 40:141, :] = (
                        +(i / time_length - 1 / 2) * tauy_trend[:, :]
                    )
                else:
                    taux_new["taux"][i, 0, 40:141, :] = (
                        +(i / time_length) * taux_trend[:, :]
                        + tauy_beg[:, :]
                    )
                    tauy_new["tauy"][i, 0, 40:141, :] = (
                        +(i / time_length) * tauy_trend[:, :]
                        + tauy_beg[:, :]
                    )

        taux_new.to_netcdf(self.setup.tau_x(it), format="NETCDF3_CLASSIC")
        tauy_new.to_netcdf(self.setup.tau_y(it), format="NETCDF3_CLASSIC")

        # doesn't change tau. TODO: Change tau
        # !!! WARNING: DOES NOTHING !!!
        # !!! JUST ADDED TO MAKE SURE PROCESS WORKS !!!
        taux_clim_obj = xr.open_dataset(self.setup.tau_clim_x(0), decode_times=False)
        taux_clim_obj.to_netcdf(
            self.setup.tau_clim_x(it),
            format="NETCDF3_CLASSIC",
        )
        tauy_clim_obj = xr.open_dataset(self.setup.tau_clim_y(0), decode_times=False)
        tauy_clim_obj.to_netcdf(
            self.setup.tau_clim_y(it),
            format="NETCDF3_CLASSIC",
        )

    def replace_dq(self, it: int) -> None:
        """
        Replace dQ variables.

        dQdf
        dQdT
        """
        dq_df_from_atm = open_dataset(self.setup.dq_output()).dq_df
        dq_df_sample = xr.open_dataarray(
            self.setup.dq_df(0),
            decode_times=False,
        )
        dq_df_new = dq_df_sample.copy()
        for t in range(12):
            dq_df_new[t, 0, 30:151, :] = can_coords(dq_df_from_atm)[:, :]
        dq_df_new.to_dataset().to_netcdf(
            self.setup.dq_df(it),
            format="NETCDF3_CLASSIC",
        )
        dq_dt_from_atm = open_dataset(self.setup.dq_output()).dq_dt
        dq_dt_sample = xr.open_dataarray(
            self.setup.dq_dt(0),
            decode_times=False,
        )
        dq_dt_new = dq_dt_sample.copy()
        for t in range(12):
            dq_dt_new[t, 0, 30:151, :] = can_coords(dq_dt_from_atm)[:, :]
        dq_dt_new.to_dataset().to_netcdf(
            self.setup.dq_dt(it),
            format="NETCDF3_CLASSIC",
        )

    def replace_surface_temp(self, it: int) -> None:
        """
        Replace sst for forcing atmosphere model.

        TODO: replace masking with actual mask.

        Args:
            it (int): iteration.
        """
        mask = open_dataset(self.setup.om_mask()).mask

        sst = can_coords(open_dataset(self.setup.om_run2f_nc()).SST_SST)
        sst_c_mean = sst.mean("T").isel(Z=0).drop("Z")

        trend_new = (
            (
                get_trend(sst + zero_Celsius, min_clim_f=True)
                .rename("ts")
                .isel(Z=0)
                .drop("Z")
            )
            .where(mask != 0.0)
            .fillna(0.0)
        )
        trend_old = xr.open_dataset(self.setup.ts_trend(0), decode_times=False)
        trend_final = trend_old.copy()
        trend_final["ts"][10:171, :] = (
            trend_final.ts[10:171, :].where(mask == 0.0).fillna(0.0)
        )
        trend_final["ts"][10:171, :] = trend_new[:, :] + trend_final.ts[10:171, :]
        # xr.testing.assert_allclose(trend_final, trend_old, atol=10)
        trend_final.fillna(0.0).to_netcdf(self.setup.ts_trend(it))

        # sst_mean: take mean
        # take mean
        sst_mean = sst_c_mean + zero_Celsius  # kelvin
        # fill in land.
        sst_mean = sst_mean.where(mask != 0.0).fillna(0.0)

        # ts_clim60
        sst_a = sst_mean.rename({"Y": "lat", "X": "lon"})
        mask_ll = mask.rename({"Y": "lat", "X": "lon"})
        sst_mean60_old = xr.open_dataset(self.setup.ts_clim60(0), decode_times=False)
        sst_mean60_final = sst_mean60_old.copy()
        sst_mean60_final["ts"][:, :] = (
            sst_a[20:141, :]
            + sst_mean60_final.ts.where(mask_ll.isel(lat=slice(20, 141)) == 0.0).fillna(
                0.0
            )[:, :]
        )
        # xr.testing.assert_allclose(sst_mean60_final, sst_mean60_old, atol=10)
        sst_mean60_final.to_netcdf(self.setup.ts_clim60(it))

        # ts_clim
        sst_b = sst_mean.rename({"Y": "Y", "X": "X"})
        sst_mean_old = xr.open_dataset(self.setup.ts_clim(0), decode_times=False)
        sst_mean_final = sst_mean_old.copy()
        sst_mean_final["ts"][10:171, :] = sst_b[:, :] + sst_mean_final.ts[
            10:171, :
        ].where(mask == 0.0).fillna(0.0)
        # xr.testing.assert_allclose(sst_mean_final, sst_mean_old, atol=10)
        sst_mean_final.to_netcdf(self.setup.ts_clim(it))

    def log(self, it: int) -> None:
        """
        Log the important information about the run.

        Args:
            it (int): Which iteration are we on?
        """
        print("logging")
        d1 = get_nino_trend(
            self.setup.om_run2f_nc(),
            self.setup.nino_png(it),
            self.setup.nino_nc(it),
        )
        d2 = get_other_trends(self.setup)
        d3 = {**d1, **d2}
        d3["it"] = it
        d3["ocean_run"] = self.ocean.run_time
        if self.cfg.wandb:
            wandb.log(d3)

    def run(self) -> None:
        """
        Run coupling.

        TODO: is this the right way to couple?
        """
        print("setting up spin up run")

        # Initial set up.
        self.ocean.compile_all()
        self.ocean.edit_run()

        if self.cfg.run:
            self.ocean.run_all(it=0)

        # atmos model.
        if self.cfg.atmos:
            # atmos takes in cfg
            self.atmos.run_all()

        self.ocean.copy_old_io(0)
        self.log(0)

        for it in range(1, self.coup.iterations):
            print(
                "coupling number ",
                it,
                " of " + str(self.coup.iterations) + " iterations.",
            )
            self.replace_dq(it)
            self.replace_stress(it)
            self.replace_surface_temp(it)
            self.ocean.edit_inputs(it)
            # self.ocean.rename(x)
            if self.cfg.run:
                self.ocean.run_all(it=it)
                self.atmos.run_all(it=it)

            # log wandb information
            self.log(it)

            # copy old io.
            self.ocean.copy_old_io(it)

        print(self.cfg.comp.sst, self.cfg.comp.prwnd, self.cfg.comp.htherm)

        plot_names = {
            "sst_"
            + str(self.cfg.comp.sst): comp_oc_sst(self.setup, str(self.cfg.comp.sst)),
            "prwnd_"
            + str(self.cfg.comp.prwnd): comp_atm_prwnd(
                self.setup, str(self.cfg.comp.prwnd)
            ),
            "htherm_"
            + str(self.cfg.comp.htherm): comp_oc_htherm(
                self.setup, str(self.cfg.comp.htherm)
            ),
        }

        # set up.s
        if self.cfg.animate:
            up_therm_qnet(self.setup, save_path=self.setup.tuq_trend_plot())
            prcp_quiver_plot(self.setup, save_path=self.setup.prcp_quiver_plot())
            self.ocean.animate_all()
            animate_coupling(self.setup)
            animate_coupling(self.setup, pac=True)
            animate_coupling(self.setup, pac=True, mask_land=True)
            animate_coupling(self.setup, pac=False, mask_land=True)
            if self.cfg.wandb:
                d_2 = {
                    "coupling_video_pac_mask_land": wandb.Video(
                        self.setup.coupling_video(pac=True, mask_land=True),
                        fps=1,
                        format="gif",
                    ),
                    "coupling_video": wandb.Video(
                        self.setup.coupling_video(pac=False, mask_land=False),
                        fps=1,
                        format="gif",
                    ),
                    "final_nino_graph": wandb.Image(
                        self.setup.nino_png(it),
                        caption=str(
                            "Final Nino region anomalies" + " over the 58 year trends"
                        ),
                    ),
                    "prcp_quiver_plot": wandb.Image(
                        self.setup.prcp_quiver_plot(),
                        caption=str(
                            "Change in precipitation and surface wind"
                            + " over the 58 years."
                        ),
                    ),
                    "tuq_trend_plot": wandb.Image(
                        self.setup.tuq_trend_plot(),
                        caption=str(
                            "Change in thermocline, upwelling and net heat flux."
                        ),
                    ),
                }

                d_3 = {}
                for i in plot_names:
                    d_3[i] = wandb.Image(plot_names[i], caption=str(i))

                if self.cfg.wandb:
                    wandb.log({**d_2, **d_3})
