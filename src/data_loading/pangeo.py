"""Get CMIP6 variables on pangeo by linking together
historical and SSP85 and historical simulations.

It is possible to change the options to link together different

Script just for monthly surface fields.

Reprocesses data onto a 1x1 degree standard lon-lat grid.
X = 0,1,2 ... 359

Rejects some models because they are difficult to process.

TODO:
   - postprocess more thorougly.

File structures inside src/data:

    Full ensemble time series:
        nc/historical.ssp585/${var}/${member}.nc

    Ensemble 60 year mean state:
        nc/historical.ssp585.mean/${var}/${member}.nc

    Ensemble 60 year climatology:
        nc/historical.ssp585.climatology/${var}/${member}.nc

    Ensemble 60 year linear trend:
        nc/historical.ssp585.trend/${var}/${member}.nc

    Ensemble 60 year nino3.4:
        nc/historical.ssp585.nino3.4/${var}.nc

    Ensemble 60 year nino3.4.climatology:
        nc/historical.ssp585.nino3.4.climatology/${var}.nc

    Ensemble 60 year nino3.4.trends:
        nc/historical.ssp585.nino3.4.trends/${var}.nc

    MMM full time series:
        nc/historical.ssp585.mmm/${var}.nc

    MMM 60 year mean:
        nc/historical.ssp585.mmm.mean/${var}.nc

    MMM 60 year trend:
        nc/historical.ssp585.mmm.trend/${var}.nc

    MMM 60 year climatology:
        nc/historical.ssp585.mmm.climatology/${var}.nc

"""
import os
from typing import Union, Callable, List, Literal, Optional
import pathlib
import numpy as np
import pandas as pd
import xarray as xr
import dask
from cmip6_preprocessing.preprocessing import combined_preprocessing
import cftime
from intake import open_catalog
import hydra
import wandb
from omegaconf import DictConfig
from src.constants import NC_PATH, CONFIG_PATH, DATA_PATH
from src.utils import timeit
from src.xr_utils import sel, can_coords, spatial_mean, get_trend, get_clim
from src.data_loading.regrid import (
    regrid_2d,
    regrid_2d_to_standard,
    regrid_1d_to_standard,
    regrid_1d,
)


# from xarray.core.variable import Variable
# This could definitely be moved to src/constants.py
START_YEAR: str = "1958"
END_YEAR: str = "2017"
DEFAULT_FUTURE_SCENARIO: str = "ssp585"

# Different scenarios that we could choose from.
SCENARIOS: List[str] = [
    "ssp119",
    "ssp126",
    "ssp245",
    "ssp370",
    "ssp434",
    "ssp460",
    "ssp585",
]

# Pangeo catalog url.
PANGEO_CAT_URL = str(
    "https://raw.githubusercontent.com/pangeo-data/"
    + "pangeo-datastore/master/intake-catalogs/master.yaml"
)

CMIP6_ENSEMBLE_VARIABLES: List[str] = [
    "pr",
    "ps",
    "clt",
    "ts",
    "sfcWind",
    "hurs",
    "tauu",
    "tauv",
]

ALL_VARIABLES: List[str] = [
    "pr",
    "ps",
    "psl",
    "clt",
    "ts",
    "sfcWind",
    "hur",
    "hurs",
    "tauu",
    "tauv",
]

# The different modelling centers.
DEFAULT_SUCCESS_LIST: List[str] = [
    "NCAR",
    "CAMS",
    "NOAA-GFDL",
    "AS-RCEC",
    "UA",
    "FIO-QLNM",
    "INM",
    "CMCC",
    "CAS",
    "CCCma",
    "NIMS-KMA",
    "NASA-GISS",
    "E3SM-Project",
    "MOHC",
    "KIOST",
    "BCC",
    "SNU",
    "CCCR-IITM",
    "THU",
]

# I think a lot of these modelling centres have problematic grids.
# Unstructured grids cannot be interpolated using xESMf,
# and both ESMpy and cdo seem to require the netcdf files to be fed into
# the algorithm, which would require a step where the full objects are
# transferred.
DEFAULT_REJECT_LIST: List[str] = [
    "AWI",
    "MRI",
    "CSIRO-ARCCSS",
    "CCCma",
    "MIROC",
    "HAMMOX-Consortium",
]

VAR_PROP_D: dict = {
    # From https://docs.google.com/spreadsheets/d/
    # 1UUtoz6Ofyjlpx5LdqhKcwHFz2SGoTQV2_yekHyMfL9Y/edit?usp=sharing
    "hurs": {
        "units": r"$\%$",
        "long_name": "Relative humidity",
        "description": "Mean relative humidity at surface",
        "limits": [0, 100],
    },
    "hur": {
        "units": r"$\%$",
        "long_name": "Relative humidity",
        "description": "Mean relative humidity at surface",
        "limits": [0, 100],
    },
    "ts": {
        "units": "K",
        "long_name": "Surface temperature",
        "description": "Surface Temperature [K]",
        "limits": [100, 373.15],  # 100K to 100 C
    },
    "vas": {
        "units": r"m s${-1}$",
        "long_name": "Northward Near-Surface Wind",
        "description": "Northward Near-Surface Wind [m s-1]",
    },
    "uas": {
        "units": r"m s${-1}$",
        "long_name": "Eastward Near-Surface Wind",
        "description": "Eastward Near-Surface Wind [m s-1]",
    },
    "pr": {
        "units": r"kg m$^{-2}$ s$^{-1}$",
        "long_name": "Precipitation",
        "description": "Precipitation [kg m-2 s-1]",
    },
    "ps": {
        "units": "Pa",
        "long_name": "Surface air pressure",
        "description": "Surface Air Pressure [Pa]",
    },
    "psl": {
        "units": "Pa",
        "long_name": "Sea level air pressure",
        "description": "Sea Level Air Pressure [Pa]",
    },
    "clt": {
        "units": r"$\%$",
        "long_name": "Total cloud cover percentage",
        "description": "Total Cloud Cover Percentage [%]",
        "limits": [0, 100],
    },
    "sfcWind": {
        "units": r"m s${-1}$",
        "long_name": "Near-Surface Wind Speed ",
        "description": "Near-Surface Wind Speed [m s-1]",
        "limits": [0, 200],  # 0 to 200 ms-1
    },
    "tauu": {
        "units": "Pa",
        "long_name": "Zonal Surface Wind Stress",
        "description": "Near-Surface Zonal Wind Stress",
    },
    "tauv": {
        "units": "Pa",
        "long_name": "Meridional Surface Wind Stress",
        "description": "Near-Surface Meridional Wind Stress",
    },
}


@np.vectorize
def standardise_time(
    time: Union[cftime.datetime, np.datetime64],
    calendar: str = "standard",  # "gregorian"
    standard_day: int = 15,
) -> cftime.datetime:
    """
    Standardise time wrapper. Makes everything start on the same day of the month.

    Args:
        time (Union[ cftime._cftime.DatetimeNoLeap, cftime._cftime.Datetime360Day,
                    np.datetime64 ]): Time array.
        calendar (str, optional): Which cftime calendar to replace it with.
            Defaults to "standard". "360_day" possible alternative.

    Returns:
        cftime._cftime.Datetime360Day: The new calendar.
    """
    if isinstance(time, np.datetime64):
        time = pd.to_datetime(time)
    # put the new time in the middle of the given month
    return cftime.datetime(
        time.year, time.month, standard_day, calendar=calendar
    )  # "360_day")


def _preproc(ds: Union[xr.Dataset, xr.DataArray]) -> Union[xr.Dataset, xr.DataArray]:
    """
    Preprocess function that includes the cmip6_preprocessing alongside custom.

    Args:
        ds (Union[xr.Dataset, xr.DataArray]): The xarray object to preprocess.

    Returns:
        Union[xr.Dataset, xr.DataArray]: The preprocessed xarray object.
    """
    dsa = ds.copy()
    dsa = combined_preprocessing(dsa)
    # dsa = dsa.assign_coords(x=(dsa.x % 360))
    dsa = dsa.assign_coords(time=("time", standardise_time(dsa.time.values)))
    return dsa


def _folder(path: str) -> None:
    """Check if folder exists, otherwise make it."""
    if not os.path.exists(path):
        os.mkdir(path)


class GetEnsemble:
    """A class to get the ensemble of CMIP6 members for monthly surface variables.

    Regrids the data on to 1 degree grid with linear 2D interpolation."""

    def __init__(
        self,
        var: str = "ts",
        table_id: str = "Amon",
        regrid: Literal["1d", "2d"] = "1d",
        output_folder: Optional[str] = None,
        regen_success_list: bool = False,
        past: str = "historical",
        future: str = DEFAULT_FUTURE_SCENARIO,
        test: bool = False,
        wandb_run: Optional[wandb.sdk.wandb_run.Run] = None,
    ) -> None:
        """
        Create the get ensemble instance and output the ensemble of netcdfs.

        TODO: add adequate unit tests to this crucial script.
        TOOD: Add example to functions.

        Args:
            var (str, optional): Variable. Defaults to "ts".
            table_id (str, optional): ID of cmip6 table to find variable in.
                Defaults to "Amon".
            regrid (Literal["1d", "2d"]):
            output_folder (Optional[str], optional): Where to output the ensemble to.
                Defaults to None, and uses default structure.
                nc/scenariomip/ts/
            regen_success_list (bool, optional): whether or not to regenerate
                the success list. Defaults to False.
            test (bool, optional): Whether or not this is a test. If it is, it
                only loads a single centre "INM". This makes tests quicker to run.
            wandb_run (Optional[wandb.sdk.wandb_run.Run]): wandb run.

        """
        self.var: str = var
        self.table_id: str = table_id
        # move to some constants file
        self.cat = open_catalog(PANGEO_CAT_URL)["climate"]["cmip6_gcs"]
        self.instit: List[str] = self.cat.unique(["institution_id"])["institution_id"][
            "values"
        ]
        self.da_lists: dict = {}
        self.regrid: Literal["1d", "2d"] = regrid
        self.past = past
        self.future = future
        if output_folder is None:
            output_folder = _folder_name(self.var, self.past + "." + self.future)
        self.output_folder = output_folder
        _folder(output_folder)
        self.wandb = wandb_run

        if test:
            self.success_list = ["INM"]
        elif regen_success_list:
            self.success_list = self.get_sucess_list()
        else:
            self.success_list = DEFAULT_SUCCESS_LIST

    @timeit
    def ensemble_timeseries(self) -> None:
        """Make a set of merged datarrays.

        Now changed to be analagous to line in Seager et al. 2019.

        Mulitimodel mean of 40 historical and RCP8.5 CMIP5 models.

        1940 to end of 2014 - historical experiment.
        2015 to end of 2099 - ssp585 experiment.

        """
        self.da_lists[self.past] = self.get_var(
            experiment=self.past,
            year_begin="1940",
            year_end="2014",
        )
        self.da_lists[self.future] = self.get_var(
            experiment=self.future, year_begin="2015", year_end="2099"
        )
        for instit in self.success_list:
            self.comp_and_match(instit=instit)

    @timeit
    def get_sucess_list(self) -> list:
        """
        The only way to generate what will work seems to be trial and error.

        Returns:
            list: success_list of model centres that I can preprocess.
        """
        failed_list = []
        empty_list = []
        success_list = []
        time_d = {}

        for i in self.instit:
            print(i)
            query = dict(
                variable_id=[self.var],
                experiment_id=[self.past],
                table_id=[self.table_id],
                institution_id=[i],
            )
            subset = self.cat.search(**query)
            z_kwargs = {"consolidated": True, "decode_times": True}

            try:
                with dask.config.set(**{"array.slicing.split_large_chunks": True}):
                    dset_dict_proc = subset.to_dataset_dict(
                        zarr_kwargs=z_kwargs, preprocess=_preproc
                    )
                if len(dset_dict_proc) == 0:
                    empty_list.append(i)
                else:
                    success_list.append(i)
                    for j in dset_dict_proc:
                        time_d[i] = dset_dict_proc[j].time.values[0]
            # pylint: disable=broad-except
            except Exception as e:
                print(e)
                print(i, "failed")
                failed_list.append(i)

        print("Success list", success_list)

        return success_list

    @timeit
    def get_var(
        self,
        experiment: str = "historical",
        year_begin: str = str(START_YEAR),
        year_end: str = "2014",
        # pylint: disable=dangerous-default-value
    ) -> xr.DataArray:
        """
        Get the variable from pangeo.

        Args:
            experiment (str, optional): Which experiment to read from.
                Defaults to "historical".
            year_begin (str, optional): The year to begin with. Defaults to "1958".
            year_end (str, optional): The year to end at. Defaults to "2014".

        Returns:
            xr.DataArray: Monthly variables with all possible ensemble members.
        """
        query = dict(
            variable_id=[self.var],
            experiment_id=[experiment],  # , "ssp585"],
            table_id=[self.table_id],
            institution_id=[
                x for x in self.success_list if x not in DEFAULT_REJECT_LIST
            ],
        )
        subset = self.cat.search(**query)

        z_kwargs = {"consolidated": True, "decode_times": True}

        # pass the preprocessing directly
        with dask.config.set(**{"array.slicing.split_large_chunks": True}):
            dset_dict_proc = subset.to_dataset_dict(
                zarr_kwargs=z_kwargs, preprocess=_preproc
            )

        da_list = []
        key_list = []

        for key in dset_dict_proc:
            da_sel = dset_dict_proc[key][self.var].sel(
                time=slice(str(year_begin), str(year_end))
            )
            if self.regrid == "2d":
                da = regrid_2d_to_standard(regrid_2d(da_sel))
            elif self.regrid == "1d":
                da = regrid_1d_to_standard(regrid_1d(da_sel))

            for i in da.member_id.values:
                key_list.append(key + "." + i)
                sub_da = da.sel(member_id=i)
                key_split = key.split(".")
                sub_da = sub_da.assign_coords(
                    {"institution": key_split[1], "model": key_split[2]}
                )
                # Try to get rid of annoying pressure / height data fields
                # from what should be a surface field.
                if "height" in sub_da.dims:
                    sub_da = sub_da.isel(height=0).drop("height")
                # pylint: disable=unnecessary-comprehension
                if "height" in [c for c in sub_da.coords]:
                    sub_da = sub_da.drop("height")
                if "plev" in sub_da.dims:
                    sub_da = sub_da.bfill("plev").isel(plev=0).drop("plev")
                # pylint: disable=unnecessary-comprehension
                if "plev" in [c for c in sub_da.coords]:
                    sub_da = sub_da.drop("plev")
                da_list.append(sub_da)

        da = xr.concat(da_list, "member")
        da = da.assign_coords({"member": key_list})  # , "time": times})

        return da

    def _coord_attrs(self, da: xr.DataArray) -> xr.DataArray:
        """Add coordinate attribute dictionaries."""
        da["X"].attrs = {
            "standard_name": "longitude",
            "units": "degrees_east",
            "long_name": "Longitude",
        }
        da["Y"].attrs = {
            "standard_name": "latitude",
            "units": "degrees_north",
            "long_name": "Latitude",
        }
        da["T"].attrs = {"standard_name": "time", "long_name": "Time"}
        return da

    def comp_and_match(self, instit: str = "NCAR") -> None:
        """
        Compare and match. Outputs the netcdf to the directory.

        Args:
            instit (str, optional): Which institute to go through. Defaults to "NCAR".
        """
        scenario_instit = self.da_lists[self.future].where(
            self.da_lists[self.future].institution == instit, drop=True
        )
        historical_instit = self.da_lists[self.past].where(
            self.da_lists[self.past].institution == instit, drop=True
        )
        # print(scenario_instit.member.values)
        # print(historical_instit.member.values)
        for model in scenario_instit.model.values:
            scenario_model = scenario_instit.where(
                scenario_instit.model == model, drop=True
            )
            historical_model = historical_instit.where(
                historical_instit.model == model, drop=True
            )
            # print(scenario_model.member.values)
            # print(historical_model.member.values)
            for member_id in scenario_instit.member_id.values:
                scenario_member = scenario_model.where(
                    scenario_model.member_id == member_id, drop=True
                )
                historical_member = historical_model.where(
                    historical_model.member_id == member_id, drop=True
                )
                if (
                    len(scenario_member.member.values) != 0
                    and len(historical_member.member.values) != 0
                ):
                    print(historical_member.member.values)
                    print(scenario_member.member.values)
                    # if len(scenario_member.member.values) !=0 and
                    # len(historical_member.member.values) != 0:
                    # if len(da_list) <= 1:
                    ensemble_da = xr.concat(
                        [
                            historical_member.isel(member=0).drop("member"),
                            scenario_member.isel(member=0).drop("member"),
                        ],
                        "T",
                    )
                    ensemble_da = self._coord_attrs(ensemble_da)
                    ensemble_da = self._wandb(ensemble_da, stage_str="processing")
                    times = ensemble_da["T"].values
                    vals, idx_start, count = np.unique(
                        times, return_counts=True, return_index=True
                    )

                    # da_list.append(ensemble_da)
                    # key_list.append(scenario_member.member.values[0])
                    ensemble_da.isel(T=idx_start).astype("float32").to_netcdf(
                        os.path.join(
                            self.output_folder,
                            self._file_name(instit, model, member_id),
                        )
                        # save some space
                        # encoding={self.var: {"dtype": "f8"}}
                    )
                else:
                    print("problem with " + model + " " + member_id)
                    print(scenario_member)
                    print(historical_member)

    def _wandb(self, da: xr.DataArray, stage_str: str = "processing") -> xr.DataArray:
        """
        ensemble_da = self._wandb(ensemble_da, stage_str = "processing")

        Args:
            da (xr.DataArray): _description_
            stage_str (str, optional): _description_. Defaults to "processing".

        Returns:
            xr.DataArray: _description_
        """
        if wandb.run is not None and self.wandb is not None:
            da.attrs["pangeo_" + stage_str + "_run"] = self.wandb.get_url()
        return da

    def _file_name(self, instit: str, model: str, member_id: str) -> str:
        return str(
            self.var
            + "."
            + _member_name(instit, model, member_id)
            + "."
            + self.past
            + "."
            + self.future
            + ".nc"
        )

    def _clip(self, da: xr.DataArray) -> xr.DataArray:
        """Clip variable between limits if they exist."""
        if "limits" in VAR_PROP_D[self.var]:
            return da.clip(
                min=VAR_PROP_D[self.var]["limits"][0],
                max=VAR_PROP_D[self.var]["limits"][1],
                keep_attrs=True,
            )
        else:
            return da

    def _mmm_attrs(self, da: xr.DataArray) -> xr.DataArray:
        """Modify variable names to make it obvious that this is a mean."""
        da.attrs["units"] = VAR_PROP_D[self.var]["units"]
        da.attrs["long_name"] = VAR_PROP_D[self.var]["long_name"] + " multi-model-mean"
        da.attrs["description"] = (
            "Multi-Model Mean " + VAR_PROP_D[self.var]["description"]
        )
        return da

    def _trend_attrs(self, da: xr.DataArray) -> xr.DataArray:
        """Modify variable names to make it obvious that this is a mean."""
        da.attrs["units"] = VAR_PROP_D[self.var]["units"]
        da.attrs["long_name"] = VAR_PROP_D[self.var]["long_name"] + " trend"
        da.attrs["description"] = "Trend for " + VAR_PROP_D[self.var]["description"]
        da.attrs["start_year"] = START_YEAR
        da.attrs["end_year"] = END_YEAR
        da = da.rename(self.var)
        return da

    def _time_mean_attrs(self, da: xr.DataArray) -> xr.DataArray:
        """Modify variable names to make it obvious that this is a mean."""
        da.attrs["units"] = VAR_PROP_D[self.var]["units"]
        da.attrs["long_name"] = VAR_PROP_D[self.var]["long_name"] + " time mean"
        da.attrs["description"] = (
            "Time mean of the " + VAR_PROP_D[self.var]["description"]
        )
        da.attrs["start_year"] = START_YEAR
        da.attrs["end_year"] = END_YEAR
        return da

    def _climatology_attrs(self, da: xr.DataArray) -> xr.DataArray:
        """Modify variable names to make it obvious that this is a mean."""
        da.attrs["units"] = VAR_PROP_D[self.var]["units"]
        da.attrs["long_name"] = VAR_PROP_D[self.var]["long_name"] + " climatology"
        da.attrs["description"] = (
            "Climatology of the " + VAR_PROP_D[self.var]["description"]
        )
        da.attrs["start_year"] = START_YEAR
        da.attrs["end_year"] = END_YEAR
        return da

    # TODO Add in minimal test before calculating MMM mean.

    def ensemble_derivatives(self) -> None:
        """Ensemble derivatives."""
        da = self.load_ensemble_da().sel(T=slice(START_YEAR, END_YEAR))
        for member in da["member"].values:
            member_da = self._clip(da.sel(member=member))
            instit = member_da.institution.values
            model = member_da.model.values
            member_id = member_da.member_id.values
            mean = self._time_mean_attrs(member_da.mean("T"))
            mean = self._wandb(mean, stage_str="time_mean")
            mean.to_netcdf(
                os.path.join(
                    _folder_name(self.var, self.past + "." + self.future + ".mean"),
                    self._file_name(instit, model, member_id),
                )
            )
            clim = self._climatology_attrs(get_clim(member_da))
            clim = self._wandb(clim, stage_str="climatology")
            clim.to_netcdf(
                os.path.join(
                    _folder_name(
                        self.var, self.past + "." + self.future + ".climatology"
                    ),
                    self._file_name(instit, model, member_id),
                )
            )
            trend = self._trend_attrs(get_trend(member_da))
            trend = self._wandb(trend, stage_str="climatology")
            trend.to_netcdf(
                os.path.join(
                    _folder_name(self.var, self.past + "." + self.future + ".trend"),
                    self._file_name(instit, model, member_id),
                )
            )

    # NINO3.4

    def _nino34_timeseries(self) -> str:
        """Multi-model mean time-series file name."""
        return os.path.join(
            _mmm_folder_name(self.past + "." + self.future + ".nino3.4"),
            self.var + ".nc",
        )

    def nino34_derivatives(self) -> None:
        """Nino3.4 derivatives."""
        da = self._clip(self.load_ensemble_da())
        minimal_members = minimal_members_ensemble()
        nino34 = spatial_mean(sel(da.sel(member=minimal_members), reg="nino3.4"))
        nino34.rename(self.var)
        nino34.to_netcdf(self._nino34_timeseries())
        clim = get_clim(nino34.sel(T=slice(START_YEAR, END_YEAR)))
        clim = clim.rename(self.var)
        clim.to_netcdf(
            os.path.join(
                _mmm_folder_name(
                    self.past + "." + self.future + ".nino3.4.climatology"
                ),
                self.var + ".nc",
            )
        )
        trends = get_60yr_trends(nino34.groupby("T.month") - clim)
        trends = trends.rename(self.var)
        trends.to_netcdf(
            os.path.join(
                _mmm_folder_name(self.past + "." + self.future + ".nino3.4.trends"),
                self.var + ".nc",
            )
        )

    def _mmm_time_series(self) -> str:
        """Multi-model mean time-series file name."""
        return os.path.join(
            _mmm_folder_name(self.past + "." + self.future + ".mmm"), self.var + ".nc"
        )

    def mmm_time_series(self):
        """Multi model mean time series."""
        da = self.load_ensemble_da()
        minimal_members = minimal_members_ensemble()
        mean = self._mmm_attrs(da.sel(member=minimal_members).mean("member"))
        mean.attrs["member_list"] = str(minimal_members)
        mean.to_netcdf(self._mmm_time_series())

    def mmm_mean_state(self) -> None:
        """
        Variable mean in time and between models for the 60 years.
        """
        mmm = xr.open_dataarray(self._mmm_time_series())
        mean = self._time_mean_attrs(mmm.sel(T=slice(START_YEAR, END_YEAR)).mean("T"))
        mean.attrs["member_list"] = mmm.attrs["member_list"]
        mean.to_netcdf(
            os.path.join(
                _mmm_folder_name(self.past + "." + self.future + ".mmm.mean"),
                self.var + ".nc",
            )
        )

    def mmm_trend(self) -> None:
        """Multi-mean trend."""
        mmm = xr.open_dataarray(self._mmm_time_series())
        trend = self._trend_attrs(get_trend(mmm.sel(T=slice(START_YEAR, END_YEAR))))
        trend.attrs["member_list"] = mmm.attrs["member_list"]
        trend.to_netcdf(
            os.path.join(
                _mmm_folder_name(self.past + "." + self.future + ".mmm.trend"),
                self.var + ".nc",
            )
        )

    def mmm_clim(self) -> None:
        """Multi-mean climatology."""
        mmm = xr.open_dataarray(self._mmm_time_series())
        clim = self._climatology_attrs(get_clim(mmm.sel(T=slice(START_YEAR, END_YEAR))))
        clim = clim.rename(self.var)
        clim.attrs["member_list"] = mmm.attrs["member_list"]
        clim.to_netcdf(
            os.path.join(
                _mmm_folder_name(self.past + "." + self.future + ".mmm.climatology"),
                self.var + ".nc",
            )
        )

    def mmm_derivatives(self) -> None:
        """Multi-model mean derivatives."""
        self.mmm_time_series()
        self.mmm_mean_state()
        self.mmm_trend()
        self.mmm_clim()

    def load_ensemble_da(self) -> xr.DataArray:
        """Load the ensemble of this variable."""
        return load_mfda(
            _folder_name(self.var, self.past + "." + self.future), self.var
        )

    def _member_list(self) -> List[str]:
        """Member list."""
        return list(self.load_ensemble_da()["member"].values)

    @timeit
    def get_future(self) -> None:
        """Get futures, make them into separate netcdfs.

        File structure:
            ScenanarioMIP/experiment_id/variable_id/ensemble_member.nc
        """
        _folder("ScenarioMIP")

        for scenario in SCENARIOS:
            _folder(os.path.join("ScenarioMIP", scenario))
            _folder(os.path.join("ScenarioMIP", scenario, self.var))
            da = self.get_var(experiment=scenario, year_begin="2014", year_end="2100")

            for member in da.member.values:

                da.sel(member=member).to_netcdf(
                    os.path.join(
                        "ScenarioMIP",
                        scenario,
                        self.var,
                        self.var
                        + "."
                        + _member_name_from_da(da.sel(member=member))
                        + ".nc",
                    )
                )


def get_60yr_trends(da: xr.DataArray) -> xr.Dataset:
    """Get 60 year trends.

    Args:
        da (xr.DataArray): input dataarray with dimension T.

    Returns:
        xr.DataArray: output dataarray with dimension end_year.
    """

    def _year(date: np.datetime64) -> int:
        year = date.astype("datetime64[Y]").astype(int) + 1970
        return int(year)

    min_year = _year(da["T"].values[0])
    max_year = _year(da["T"].values[-1])
    beginning = min_year + 59
    finish = max_year
    da_tr_list = []
    if beginning < finish:
        for start, end in [(str(x - 59), str(x)) for x in range(beginning, finish + 1)]:
            da_tr_list.append(
                get_trend(da.sel(T=slice(start, end)), keep_ds=True).expand_dims(
                    {"end_year": [int(end)]}
                )
            )
    da_tr = xr.merge(da_tr_list)
    da_tr["end_year"].attrs["long_name"] = "End year"
    da_tr["rise"].attrs["units"] = "K"
    da_tr["rise"].attrs["long_name"] = "Trend over 60 years up to end year"
    return da_tr.rise


def load_ensemble_da(var: str, path: str) -> xr.DataArray:
    """
    Args:
        var (str): the variable to load.
            E.g. "ts"
        path (str): the path to the variable to load.
            E.g. "src/data/nc/nc_ts"

    Returns:
        xr.DataArray: the xarray datarray with "X", "Y", "T", and "member"
            as dimensions.
    """
    return load_mfda(path, var)


def _scenariomip_data(var: str = "ts", scenario: str = "ssp585"):
    return os.path.join(os.path.join("ScenarioMIP", scenario, var))


def _print_scenario_folders():
    for scenario in SCENARIOS:
        print(_scenariomip_data(scenario=scenario))


def _scenario_nino34(scenario: str = "ssp585"):
    return load_mfda(_scenariomip_data(scenario=scenario, var="ts_nino3.4"), "ts")


def _scenario_da(var: str = "ts", scenario: str = "ssp585"):
    return load_mfda(_scenariomip_data(scenario=scenario, var=var), "ts")


def load_mfds(files_path: str, var: str) -> xr.Dataset:
    """
    Load multi file dataset.

    Args:
        files_path (str): Path to files.
        var (str): variable.

    Returns:
        xr.Dataset: The variable dataset.
    """
    return xr.open_mfdataset(
        os.path.join(files_path, "*.nc"),
        concat_dim="member",
        preprocess=_member_preproc_func(var),
    )


def load_mfda(files_path: str, var: str) -> xr.DataArray:
    """
    Load multi file dataarray.

    Args:
        files_path (str): Path to files.
        var (str): variable.

    Returns:
        xr.DataArray: The variable dataset.
    """
    return load_mfds(files_path, var)[var]


def _print_scenario_sel():
    for scenario in SCENARIOS:
        print(can_coords(load_mfda(_scenariomip_data(scenario=scenario), "ts")))


def _member_name_from_da(da: xr.DataArray) -> str:
    """Member name from da already referenced by member.

    Args:
        da (xr.DataArray): xarray input.

    Returns:
        str: member name.
    """
    instit = da.institution.values
    model = da.model.values
    member_id = da.member_id.values
    return _member_name(instit, model, member_id)


def _member_name(instit: str, model: str, member_id: str) -> str:
    """
    Member name for dimension, and for file names.

    Args:
        instit (str): CMIP6 modelling institution.
        model (str): Model name.
        member_id (str): Model run id.

    Returns:
        str: ${istit}.${model}.${member_id}
    """
    return str(instit) + "." + str(model) + "." + str(member_id)


def _folder_name(
    var: str, experiment: str, root_direc: Union[pathlib.Path, str] = NC_PATH
) -> str:
    """
    Folder name for ensemble of netcdfs.

    Args:
        var (str): Variable name (e.g "pr").
        experiment (str): Experiment name (e.g. "historical")
        root_direc (Union[pathlib.Path, str]): Folder path.

    Returns:
        str: Folder path (e.g. "pr")
    """
    _folder(root_direc)
    experiment_direc = str(os.path.join(root_direc, experiment))
    _folder(experiment_direc)
    path = str(os.path.join(experiment_direc, var))
    _folder(path)
    return path


def _mmm_folder_name(
    experiment: str, root_direc: Union[pathlib.Path, str] = NC_PATH
) -> str:
    """
    Folder name for ensemble of netcdfs.

    Args:
        experiment (str): Experiment name (e.g. "historical")
        root_direc (Union[pathlib.Path, str]): Folder path.

    Returns:
        str: Folder path (e.g. "pr")
    """
    _folder(root_direc)
    experiment_direc = str(os.path.join(root_direc, experiment))
    _folder(experiment_direc)
    return experiment_direc


def _member_preproc_func(var_str: str) -> Callable:
    """Preprocessing function for loading preprocessed data.

    Adds the member dimension and coordinate.
    """

    def _preproc_func(ds: xr.Dataset) -> xr.Dataset:
        dsa = ds.copy()
        dsa = dsa.expand_dims("member")
        if "__xarray_dataarray_variable__" in dsa:
            dsa = dsa["__xarray_dataarray_variable__"].rename(var_str).to_dataset()
        dsa = dsa.assign_coords(
            {
                "member": [
                    _member_name(
                        str(dsa.institution.values),
                        str(dsa.model.values),
                        str(dsa.member_id.values),
                    )
                ]
            }
        )
        return dsa

    return _preproc_func


def make_future_nino34(
    scenario: Literal["ssp126", "ssp245", "ssp370", "ssp585"] = "ssp585"
) -> None:
    """
    Make future `nino3.4` records.

    Args:
        scenario (Literal["ssp126", "ssp245", "ssp370", "ssp585"]):
            Which scenario to get data from.
    """
    da = _scenario_da(scenario=scenario)
    new_da = regrid_2d_to_standard(da)
    nino34 = spatial_mean(sel(new_da, reg="nino3.4"))
    output_folder = _scenariomip_data(var="ts_nino3.4", scenario=scenario)
    _folder(output_folder)
    for i in range(len(nino34.member.values)):
        da_i = nino34.isel(member=i)
        member_name = _member_name_from_da(da_i)
        output_file_name = os.path.join(
            output_folder, "ts_nino3.4." + member_name + ".nc"
        )
        da_i.to_netcdf(output_file_name)


def plausible_temperatures(values: np.ndarray) -> bool:
    """
    Checks if plausible temperatures in Kelvin.

    Args:
        values (np.ndarray): numpy array.

    Returns:
        bool: Whether or not the temperatures are all plausible.
    """
    lower_bound = 100  # 100K
    upper_bound = 373.15  # 100 Degrees Celsius
    too_high = values > lower_bound
    too_low = values < upper_bound
    not_nan = np.invert(np.isnan(values))
    return np.all(too_high) and np.all(too_low) and np.all(not_nan)


def test_if_regridding_ok() -> None:
    """Trying to test if different regridding methods each work."""
    import matplotlib.pyplot as plt

    oned = GetEnsemble(regrid="1d", test=True)
    twod = GetEnsemble(regrid="2d", test=True)
    da_list = [oned.get_var(), twod.get_var()]
    for i, da in enumerate(da_list):
        print(da)
        print(da.isel(T=0, member=0).values)
        print(plausible_temperatures(da.isel(T=0, member=0).values))
        assert not plausible_temperatures(np.array([0, 200]))
        da.isel(T=0, member=0).plot()
        plt.savefig("fig-" + str(i) + "-a.png")
        plt.clf()


@timeit
def members_df() -> pd.DataFrame:
    """
    Members dataframe showing whether each variable exists.

    Returns:
        pd.DataFrame: Truth/False for each variable for each possible ensemble member.
    """
    var_dict = {}
    truth_dict = {}
    variables: List[str] = ALL_VARIABLES
    combined_list = []
    for var in variables:
        ens = GetEnsemble(
            var=var,
            regen_success_list=False,
        )
        var_dict[var] = ens._member_list()
        combined_list = [*combined_list, *var_dict[var]]

    unique_list = list(np.unique(combined_list))

    for var in var_dict:
        truth_dict[var] = [x in var_dict[var] for x in unique_list]

    return pd.DataFrame(data=truth_dict, index=unique_list)


def minimal_ensemble(df: pd.DataFrame) -> List[str]:
    """
    Minimal ensemble.

    Args:
        df (pd.DataFrame): members_df.

    Returns:
        List[str]: list of ensemble members
    """
    df = df[CMIP6_ENSEMBLE_VARIABLES]
    print(df)
    output = []
    for index, row in df.iterrows():
        if np.all(row):
            output.append(index)
    print(len(output))
    return output


def minimal_members_ensemble() -> List[str]:
    """
    Minimal members ensemble.

    Returns:
        List[str]: Ensemble member name string list.
    """
    return minimal_ensemble(
        pd.read_csv(DATA_PATH / "ensemble_variable_members.csv", index_col=0)
    )


def members_csv() -> None:
    """Make members csv so that I can transfer it."""
    members_df().to_csv(DATA_PATH / "ensemble_variable_members.csv")


@timeit
@hydra.main(config_path=CONFIG_PATH, config_name="pangeo.yml")
def main(cfg: DictConfig) -> None:
    """
    Run the data generation scripts.

    Args:
        cfg (DictConfig): The config struct to run the program.
    """
    print(cfg)
    wandb_run = wandb.init(
        project=cfg.project,
        entity=cfg.user,
        save_code=True,
        name=cfg.name,
        # pylint: disable=protected-access
        config=cfg._content,
    )
    ens = GetEnsemble(
        var=cfg.var,
        regen_success_list=cfg.regen_success_list,
        wandb_run=wandb_run,
    )
    if cfg.ensemble_timeseries:
        ens.ensemble_timeseries()
    if cfg.regenerate_members_csv:
        members_csv()
    if cfg.ensemble_derivatives:
        print(ens._member_list())
        ens.ensemble_derivatives()
    if cfg.mmm_derivatives:
        ens.mmm_derivatives()
    if cfg.nino34_derivatives:
        ens.nino34_derivatives()
    wandb.finish()


if __name__ == "__main__":
    # from src.data_loading.pangeo import GetEnsemble
    # python src/data_loading/pangeo.py > clt.txt
    # for var_str in VAR_PROP_D:
    #    GetEnsemble(
    #        var=var_str, output_folder=_folder_name(var_str), regen_success_list=True
    #    )
    # minimal_members_ensemble()
    # for var_str in VAR_PROP_D:
    #    mean_var(var=var_str)
    # make_wsp()
    # from src.data_loading.pangeo import get_vars
    # get_vars(["hurs", "psl"], )
    # get_vars(["clt"])
    # ts = GetEnsemble(var="ts")
    # ts.get_future()
    # _print_scenario_sel()
    # for scenario in SCENARIOS:
    #    make_future_nino34(scenario=scenario)
    # test_if_regridding_ok()
    # get_vars(["ts", "sfcWind", "hurs"], regen_success_list=False)
    # pylint: disable=no-value-for-parameter
    main()
    # print(members_df())
    # members_csv()
    # conda activate pangeo3
    # python src/data_loading/pangeo.py -m var=sfcWind,hurs
    # python src/data_loading/pangeo.py -m var=ps,ts
    # python src/data_loading/pangeo.py -m var=clt,pr
    # python src/data_loading/pangeo.py -m var=tauv,tauu
    # python src/data_loading/pangeo.py -m var=pr,ps,clt,ts,sfcWind,hurs ensemble_timeseries=false
    # python src/data_loading/pangeo.py -m var=pr,ps,clt,ts,sfcWind,hurs ensemble_timeseries=false ensemble_derivatives=false mmm_derivatives=false
