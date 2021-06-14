"""Coupling between ocean and atmos models is through the sea surface stress.

Example:
    Import statement usage::
        from src.models.coupling import Coupling

"""
import os
import shutil
from typing import Tuple
import xarray as xr
from typeguard import typechecked
from omegaconf import DictConfig
from src.models.model_setup import ModelSetup
from src.models.atmos import Atmos
from src.models.ocean import Ocean
from src.xr_utils import can_coords, open_dataset, cut_and_taper, get_trend


# pylint: disable=no-value-for-parameter
class Coupling:
    """
    Coupled model part.

    Model solution method. The atmosphere equations are solved by Fourier
    transforming in longitude, forming an equation for v for each zonal wavenumber
    that is finite differenced, and the resulting tri-diagonal system
    is solved by matrix
    inversion, transforming back into longitude. Finally, u and Φ are derived by
    backsubstitution. The ocean equations are solved using the ‘INC’ scheme31,
    integrating the model forward, after spin-up with climatological conditions,
    forced by the time-varying ECMWF wind stress and, for the case with CO2 forcing,
    changing f′1
    in the net surface longwave radiation calculation. Change over 1958–2017
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
    value of CO2. The model wind stress change is computed as ρ c Wu a D , where cD
    is a drag coefficient and u is the vector surface wind change computed by the
    atmosphere model, which is added to the ECMWF climatological stresses. Since
    the atmosphere model dynamics are only applicable in the tropics, the computed
    wind stress anomaly is only applied to the ocean model between 20 S and 20 N,
    and is linearly tapered to zero at 25 S and 2 5N.
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
        wind_speed_mean: float,
        u_wind: xr.DataArray,
        v_wind: xr.DataArray,
    ) -> Tuple[xr.DataArray, xr.DataArray]:
        """Wind stress flux.

        Args:
            wind_speed_mean (float): the climatological annual
                mean wind speed, which is taken from ECMWF
                reanalysis for our standard model and from the CMIP5
                multimodel mean when examining causes of bias in the
                CMIP5 model.
            u_wind (xr.DataArray): The zonal wind speed.
            v_wind (xr.DataArray): The meridional wind speed.

        Returns:
            Tuple[xr.DataArray, xr.DataArray]: [zonal wind stress,
                meridional wind stress]

        """
        stress_coeff = self.coup.rho_air * self.coup.c_d * wind_speed_mean
        return stress_coeff * u_wind, stress_coeff * v_wind

    def get_wind_speed_mean(self, file_name: str = "") -> float:
        """Get wind speed mean."""
        print("get wind speed mean")
        xr.open_dataset(file_name).mean("T")

    def replace_dq(self, it: int) -> None:
        """
        Replace dQ variables.

        dQdf
        dQdT
        """
        dq_df_from_atm = open_dataset(self.setup.dq_output()).dq_df
        dq_df_sample = xr.open_dataarray(
            os.path.join(self.setup.ocean_data_path, "dQdf-sample.nc"),
            decode_times=False,
        )
        dq_df_new = dq_df_sample.copy()
        for t in range(12):
            dq_df_new[t, 0, 30:151, :] = can_coords(dq_df_from_atm)[:, :]
        dq_df_new.to_dataset().to_netcdf(
            self.setup.dq_df(it),
            format="NETCDF3_CLASSIC",
        )
        dq_dt_from_atm = open_dataset(
            os.path.join(self.setup.atmos_path, "dQ.nc")
        ).dq_dt
        dq_dt_sample = xr.open_dataarray(
            os.path.join(self.setup.ocean_data_path, "dQdT-sample.nc"),
            decode_times=False,
        )
        dq_dt_new = dq_dt_sample.copy()
        for t in range(12):
            dq_dt_new[t, 0, 30:151, :] = can_coords(dq_dt_from_atm)[:, :]
        dq_dt_new.to_dataset().to_netcdf(
            self.setup.dq_dt(it),
            format="NETCDF3_CLASSIC",
        )

    def replace_stress(self, it: int) -> None:
        """Replace the stress files. Currently just resaves the files."""

        taux_obj = xr.open_dataset(self.setup.tau_x(0), decode_times=False)
        taux_obj.to_netcdf(
            self.setup.tau_x(it),
            format="NETCDF3_CLASSIC",
        )
        tauy_obj = xr.open_dataset(self.setup.tau_y(0), decode_times=False)
        tauy_obj.to_netcdf(
            self.setup.tau_y(it),
            format="NETCDF3_CLASSIC",
        )

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

        # test the cut and taper functions work.
        cut_and_taper(can_coords(xr.open_dataset(self.setup.tcam_output()).vtrend))
        cut_and_taper(can_coords(xr.open_dataset(self.setup.tcam_output()).utrend))

    def replace_surface_temp(self, it: int) -> None:
        """
        Replace sst for forcing atmosphere model.

        Args:
            it (int): iteration.
        """

        print(it)
        print(self.setup.ts_trend(it - 1))
        print(self.setup.ts_clim(it - 1))
        print(self.setup.ts_clim60(it - 1))
        print(self.setup.ts_trend(it))
        print(self.setup.ts_clim(it))
        print(self.setup.ts_clim60(it))
        shutil.copy(self.setup.ts_trend(0), self.setup.ts_trend(it))
        shutil.copy(self.setup.ts_clim(0), self.setup.ts_clim(it))
        shutil.copy(self.setup.ts_clim60(0), self.setup.ts_clim60(it))

        sst = can_coords(open_dataset(self.setup.om_run2f_nc()).SST_SST)
        # sst = sst.where(sst != 0.0)
        trend_new = (
            get_trend(sst + self.cfg.atm.temp_0_c, min_clim_f=True)
            .rename("ts")
            .isel(Z=0)
            .drop("Z")
        )
        trend_old = xr.open_dataset(self.setup.ts_trend(0), decode_times=False).ts
        trend_final = trend_old.copy()
        trend_final[10:171, :] = trend_new[:, :]
        trend_fin_ds = trend_final.to_dataset(name="ts")
        trend_fin_ds.to_netcdf(self.setup.ts_trend(it))

    def run(self) -> None:
        """
        Run coupling.
        """
        print("setting up spin up run")

        # Initial set up.
        self.ocean.compile_all()
        self.ocean.edit_run()

        print("run")
        if self.cfg.run:
            self.ocean.run_all(it=0)

        # atmos model.
        if self.cfg.atmos:
            # atmos takes in cfg
            self.atmos.run_all()

        self.ocean.copy_old_io(0)

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
                self.ocean.copy_old_io(it)

        # set up.
        if self.cfg.animate:
            self.ocean.animate_all()
