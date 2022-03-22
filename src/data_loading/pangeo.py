"""Get CMIP6 variables on pangeo by linking together
historical and SSP85 and historical simulations.

Script just for monthly surface fields.

Reprocesses data onto a 1x1 degree standard lon-lat grid.

Rejects some models because they are difficult to process.

TODO:
   - Grid keyword
   - postprocess
"""
import os
from typing import Union, Callable, List, Literal
import numpy as np
import pandas as pd
import xarray as xr
import xesmf as xe
import dask
from cmip6_preprocessing.preprocessing import combined_preprocessing
import cftime
from intake import open_catalog
from src.utils import timeit
from src.xr_utils import sel, can_coords, spatial_mean
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
FUTURE_SCENARIO = "ssp585"

SCENARIOS: List[str] = ["ssp585", "ssp245", "ssp370", "ssp126"]

PANGEO_CAT_URL = str(
    "https://raw.githubusercontent.com/pangeo-data/"
    + "pangeo-datastore/master/intake-catalogs/master.yaml"
)

DEFAULT_SUCCESS_LIST = [
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
DEFAULT_REJECT_LIST = [
    "AWI",
    "MRI",
    "CSIRO-ARCCSS",
    "CCCma",
    "MIROC",
    "HAMMOX-Consortium",
]

VAR_PROP_D = {
    # From https://docs.google.com/spreadsheets/d/
    # 1UUtoz6Ofyjlpx5LdqhKcwHFz2SGoTQV2_yekHyMfL9Y/edit?usp=sharing
    # psl - sea level pressure.
    # hurs - Near-Surface Relative Humidity
    "hurs": {
        "units": r"$\%$",
        "long_name": "Relative humidity",
        "description": "Mean relative humidity at surface",
    },
    "hur": {
        "units": r"$\%$",
        "long_name": "Relative humidity",
        "description": "Mean relative humidity at surface",
    },
    "ts": {
        "units": "K",
        "long_name": "Suface temperature",
        "description": "Surface Temperature [K]",
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
        "long_name": "Suface air pressure",
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
    },
    "sfcWind": {
        "units": r"m s${-1}$",
        "long_name": "Near-Surface Wind Speed ",
        "description": "Near-Surface Wind Speed [m s-1]",
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
    Preprocess.

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
        output_folder: str = "nc",
        regen_success_list: bool = False,
        test: bool = False,
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
            output_folder (str, optional): Where to output the ensemble to.
                Defaults to "nc".
            regen_success_list (bool, optional): whether or not to regenerate
                the success list. Defaults to False.
            test (bool, optional): Whether or not this is a test. If it is, it
                only loads a single centre "INM".

        """
        self.output_folder = output_folder
        self.var = var
        self.table_id = table_id
        # move to some constants file
        self.cat = open_catalog(PANGEO_CAT_URL)["climate"]["cmip6_gcs"]
        self.instit: List[str] = self.cat.unique(["institution_id"])["institution_id"][
            "values"
        ]
        self.da_lists: dict = {}
        self.regrid: Literal["1d", "2d"] = regrid

        if test:
            self.success_list = ["INM"]
        elif regen_success_list:
            self.success_list = self.get_sucess_list()
        else:
            self.success_list = DEFAULT_SUCCESS_LIST

    @timeit
    def make_comb80_das(self) -> None:
        """Make a set of merged datarrays."""
        self.da_lists["hist"] = self.get_var(
            experiment="historical",
            year_begin="1948",
            year_end="2018",
        )
        self.da_lists["ssp585"] = self.get_var(
            experiment="ssp585", year_begin="2014", year_end="2027"
        )
        for instit in self.success_list:
            self.comp_and_match(instit=instit)
        # comp_and_match(instit="NCAR")
        # print(da_list)

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
                experiment_id=["historical"],  # , "ssp585"],
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
                print("regridding 2d")
                da = regrid_2d_to_standard(regrid_2d(da_sel))
            elif self.regrid == "1d":
                print("regridding 1d")
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

    def comp_and_match(self, instit: str = "NCAR") -> None:
        """
        Compare and match. Outputs the netcdf to the directory.

        Args:
            instit (str, optional): Which institute to go through. Defaults to "NCAR".
        """
        ssp585_instit = self.da_lists["ssp585"].where(
            sself.da_lists["ssp585"].institution == instit, drop=True
        )
        hist_instit = self.da_lists["hist"].where(
            self.da_lists["hist"].institution == instit, drop=True
        )
        # print(ssp585_instit.member.values)
        # print(hist_instit.member.values)
        for model in ssp585_instit.model.values:
            ssp585_model = ssp585_instit.where(ssp585_instit.model == model, drop=True)
            hist_model = hist_instit.where(hist_instit.model == model, drop=True)
            # print(ssp585_model.member.values)
            # print(hist_model.member.values)
            for member_id in ssp585_instit.member_id.values:
                ssp585_member = ssp585_model.where(
                    ssp585_model.member_id == member_id, drop=True
                )
                hist_member = hist_model.where(
                    hist_model.member_id == member_id, drop=True
                )
                if (
                    len(ssp585_member.member.values) != 0
                    and len(hist_member.member.values) != 0
                ):
                    print(hist_member.member.values)
                    print(ssp585_member.member.values)
                    # if len(ssp585_member.member.values) !=0 and
                    # len(hist_member.member.values) != 0:
                    # if len(da_list) <= 1:
                    ensemble_da = xr.concat(
                        [
                            hist_member.isel(member=0).drop("member"),
                            ssp585_member.isel(member=0).drop("member"),
                        ],
                        "time",
                    )
                    # da_list.append(ensemble_da)
                    # key_list.append(ssp585_member.member.values[0])
                    ensemble_da.to_netcdf(
                        os.path.join(
                            self.output_folder,
                            self.var
                            + "."
                            + _member_name(instit, model, member_id)
                            + ".80.nc",
                        )
                    )
                else:
                    print("problem with " + model + " " + member_id)
                    print(ssp585_member)
                    print(hist_member)
        # da = xr.concat(da_list, "member")
        # da = da.assign_coords({"member": key_list})
        # return da_list

    def mean_var(self) -> None:
        """
        Variable mean in time and between models for the 60 years.

        Args:
            var (str, optional): Variable name string. Defaults to "ts".
        """
        da = xr.open_mfdataset(
            _folder_name(self.var) + "/*.nc",
            concat_dim="member",
            preprocess=_get_preproc_func(var),
        )
        print(da)
        mean = da.sel(time=slice(START_YEAR, END_YEAR)).mean("member").mean("time")
        mean[var].attrs["units"] = VAR_PROP_D[var]["units"]
        mean[var].attrs["long_name"] = VAR_PROP_D[var]["long_name"] + " mean"
        mean[var].attrs["description"] = "Mean " + VAR_PROP_D[var]["description"]
        _folder(_folder_name("mean"))
        mean.to_netcdf(os.path.join(_folder_name("mean"), var + ".nc"))

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


def _scenariomip_data(var: str = "ts", scenario: str = "ssp585"):
    return os.path.join(os.path.join("ScenarioMIP", scenario, var))


def _print_scenario_folders():
    for scenario in SCENARIOS:
        print(_scenariomip_data(scenario=scenario))


def _scenario_nino34(scenario: str = "ssp585"):
    return xr.open_mfdataset(
        os.path.join(_scenariomip_data(scenario=scenario, var="ts_nino3.4"), "*.nc"),
        concat_dim="member",
        preprocess=_get_preproc_func("ts"),
    )["ts"]


def _scenario_da(var: str = "ts", scenario: str = "ssp585"):
    return xr.open_mfdataset(
        os.path.join(_scenariomip_data(scenario=scenario, var=var), "*.nc"),
        concat_dim="member",
        preprocess=_get_preproc_func(var),
    )[var]


def _print_scenario_sel():
    for scenario in SCENARIOS:
        print(
            can_coords(
                xr.open_mfdataset(
                    os.path.join(_scenariomip_data(scenario=scenario), "*.nc"),
                    concat_dim="member",
                    preprocess=_get_preproc_func("ts"),
                ).ts
            )
        )


def _member_name_from_da(da: xr.DataArray) -> str:
    """member name from da already referenced by member"""
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
        str: thing.
    """
    return str(instit) + "." + str(model) + "." + str(member_id)


def _folder_name(var: str) -> str:
    """
    Folder name for ensemble of netcdfs.

    Args:
        var (str): Variable name (e.g "pr")

    Returns:
        str: Folder path (e.g. "nc_pr")
    """
    return "nc_" + var


def _get_preproc_func(var_str: str) -> Callable:
    """Preprocessing function for later 'make_${resource}' functions."""

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


def get_vars(var_list: List[str], regen_success_list=False) -> None:
    """
    Get the variable means and ensembles for the enviroment.

    Args:
        var_list (List[str]): List of variables to make.
    """
    for var_str in var_list:
        ens = GetEnsemble(
            var=var_str,
            output_folder=_folder_name(var_str),
            regen_success_list=regen_success_list,
        )
        ens.make_das()
        ens.mean_var()


def make_future_nino34(scenario="ssp585") -> None:
    da = _scenario_da(scenario=scenario)
    new_da = regridded_to_standard(da)
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
    Cheks if plausible temperatures in Kelvin.

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
    """Trying to test if different regridding mehtods each work."""
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
        plt.savefig("fig" + str(i) + ".png")
        plt.clf()


if __name__ == "__main__":
    # from src.data_loading.pangeo import GetEnsemble
    # python src/data_loading/pangeo.py > clt.txt
    # for var_str in VAR_PROP_D:
    #    GetEnsemble(
    #        var=var_str, output_folder=_folder_name(var_str), regen_success_list=True
    #    )
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
    test_if_regridding_ok()
