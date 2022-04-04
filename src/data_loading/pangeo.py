"""Get CMIP6 variables on pangeo by linking together
historical and SSP85 and historical simulations."""
import os
from typing import Union, Callable, List
import numpy as np
import pandas as pd
import xarray as xr
import dask
from cmip6_preprocessing.preprocessing import combined_preprocessing
import cftime
from intake import open_catalog

# from xarray.core.variable import Variable

START_YEAR: str = "1958"
END_YEAR: str = "2017"

DEFAULT_SUCCESS_LIST = [
    # "NCAR",
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
    # "CCCR-IITM",
    "THU",
]

PANGEO_CAT_URL = str(
    "https://raw.githubusercontent.com/pangeo-data/"
    + "pangeo-datastore/master/intake-catalogs/master.yaml"
)

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
    Standardise time.

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
    dsa = dsa.assign_coords(time=standardise_time(dsa.time.values))
    return dsa


def _folder(path: str) -> None:
    """Check if folder exists, otherwise make it."""
    if not os.path.exists(path):
        os.mkdir(path)


class GetEnsemble:
    """A class to get the ensemble of CMIP6 members for monthly surface variables.

    Regrids the data on to 1 degree grid with linear 2D interpolation."""

    def __init__(
        self, var: str = "ts", output_folder: str = "nc", regen_success_list=False
    ) -> None:
        """
        Create the get ensemble instance and output the ensemble of netcdfs.

        Args:
            var (str, optional): Variable. Defaults to "ts".
            output_folder (str, optional): Where to output the ensemble to.
                Defaults to "nc".
            regen_success_list (bool, optional): whether or not to regenerate
                the success list. Defaults to False.

        """
        self.output_folder = output_folder
        _folder(output_folder)

        self.var = var
        # move to some constants file
        self.cat = open_catalog(PANGEO_CAT_URL)["climate"]["cmip6_gcs"]
        self.instit = self.cat.unique(["institution_id"])["institution_id"]["values"]

        if regen_success_list:
            self.success_list = self.get_sucess_list()
        else:
            self.success_list = DEFAULT_SUCCESS_LIST
        self.da_hist = self.get_var(
            xlim=[0, 360],
            ylim=[-80, 80],
            var=self.var,
            year_begin="1948",
        )
        self.da_ssp585 = self.get_var(
            experiment="ssp585",
            year_begin="2014",
            year_end="2027",
            xlim=[0, 360],
            ylim=[-80, 80],
            var=self.var,
        )
        for instit in self.success_list:
            self.comp_and_match(instit=instit)
        # comp_and_match(instit="NCAR")
        # print(da_list)

    def change_t_axis(
        self,
        ds: xr.Dataset,
        calendar: str = "standard",  # "gregorian"
    ) -> xr.Dataset:
        """
        Change the time axis.

        Args:
            ds (xr.Dataset): The dataset to change.
            calendar (str, optional): [description]. Defaults to "standard".

        Returns:
            xr.Dataset: The dataset with the new time axis.
        """
        # ds = change_t_axis(ds, calendar="360_day")
        dsa = ds.copy()
        return dsa.assign_coords(
            time=standardise_time(ds.time.values), calendar=calendar
        )

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
                table_id=["Amon"],
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

        return success_list

    def get_var(
        self,
        experiment: str = "historical",
        year_begin: str = START_YEAR,
        year_end: str = "2014",
        var: str = "ts",
        # pylint: disable=dangerous-default-value
        xlim: List[int] = [100, 290],
        ylim: List[int] = [-30, 30],
    ) -> xr.DataArray:
        """
        Get the variable from pangeo.

        Args:
            experiment (str, optional): Which experiment to read from.
                Defaults to "historical".
            year_begin (str, optional): The year to begin with. Defaults to "1958".
            year_end (str, optional): The year to end at. Defaults to "2014".
            var (str, optional): The variable to select. Defaults to "ts".
            xlim (list, optional): The longitude limits. Defaults to [100, 290].
            ylim (list, optional): The latitude limits. Defaults to [-30, 30].

        Returns:
            xr.DataArray: Monthly variables with all possible ensemble members.
        """
        query = dict(
            variable_id=[var],
            experiment_id=[experiment],  # , "ssp585"],
            table_id=["Amon"],
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
            print(key)
            da = (
                dset_dict_proc[key][var]
                .sel(
                    x=slice(xlim[0] - 1, xlim[1] + 1),
                    y=slice(ylim[0] - 1, ylim[1] + 1),
                    time=slice(year_begin, year_end),
                )
                .interp(
                    x=list(range(xlim[0], xlim[1] + 1)),
                    y=list(range(ylim[0], ylim[1] + 1)),
                )
            )

            for i in da.member_id.values:
                key_list.append(key + "." + i)
                sub_da = da.sel(member_id=i)
                key_split = key.split(".")
                sub_da = sub_da.assign_coords(
                    {"institution": key_split[1], "model": key_split[2]}
                )
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

    def get_ensemble(self, var: str = "ts") -> xr.DataArray:
        """
        Get a variable ensemble.

        Args:
            var (str, optional): The variable. Defaults to "ts".

        Returns:
            xr.DataArray: Ensemble datarray full of different items.
        """
        da_hist = self.get_var(xlim=[0, 360], ylim=[-80, 80], var=var)
        da_ssp585 = self.get_var(
            experiment="ssp585",
            year_begin="2014",
            year_end=END_YEAR,
            xlim=[0, 360],
            ylim=[-80, 80],
            var=var,
        )
        da_comb = xr.concat(
            [da_hist.mean("model_center"), da_ssp585.mean("model_center")], "time"
        )
        da_comb.attrs["hist_list"] = str(da_hist["model_center"].values)
        da_comb.attrs["ssp585_list"] = str(da_ssp585["model_center"].values)
        return da_comb

    def comp_and_match(self, instit: str = "NCAR") -> None:
        """
        Compare and match.

        Args:
            instit (str, optional): Which institute to go through. Defaults to "NCAR".
        """
        ssp585_instit = self.da_ssp585.where(
            self.da_ssp585.institution == instit, drop=True
        )
        hist_instit = self.da_hist.where(self.da_hist.institution == instit, drop=True)
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
        # da = xr.concat(da_list, "member")
        # da = da.assign_coords({"member": key_list})
        # return da_list


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


def mean_var(var: str = "ts") -> None:
    """
    Variable mean in time and between models.

    Args:
        var (str, optional): Variable name string. Defaults to "ts".
    """
    da = (
        xr.open_mfdataset(
            _folder_name(var) + "/*.nc",
            concat_dim="member",
            preprocess=_get_preproc_func(var),
        )
        .drop("lon")
        .drop("lat")
    )
    print(da)
    mean = da.sel(time=slice(START_YEAR, END_YEAR)).mean("member").mean("time")
    mean[var].attrs["units"] = VAR_PROP_D[var]["units"]
    mean[var].attrs["long_name"] = VAR_PROP_D[var]["long name"] + " mean"
    mean[var].attrs["description"] = "Mean " + VAR_PROP_D[var]["description"]
    mean.to_netcdf(os.path.join(_folder_name("mean"), var + ".nc"))


def make_wsp() -> None:
    """
    Make wsp from uas and vas.
    """
    _folder("nc_wsp")

    for u, v in [
        (
            sorted(os.listdir(_folder_name("uas")))[i],
            sorted(os.listdir(_folder_name("vas")))[i],
        )
        for i in range(len(os.listdir(_folder_name("uas"))))
    ]:
        print(
            os.path.join(_folder_name("uas"), u), os.path.join(_folder_name("vas"), v)
        )
        wsp = np.sqrt(
            np.square(
                xr.open_dataarray(os.path.join(_folder_name("uas"), u))
                .drop("lon")
                .drop("lat")
            )
            + np.square(
                xr.open_dataarray(os.path.join(_folder_name("vas"), v))
                .drop("lon")
                .drop("lat")
            )
        ).ffill("x")
        wsp.attrs["long_name"] = "Mean wind speed at near surface"
        w = "wsp." + u.strip("uas.")
        wsp.rename("wsp").to_netcdf(os.path.join(_folder_name("wsp"), w))


def mean_wsp() -> None:
    """
    Make mean windspeed.
    """
    wsp_da = xr.open_mfdataset(
        _folder_name("wsp") + "/*.nc",
        concat_dim="member",
        preprocess=_get_preproc_func("wsp"),
    )
    wsp_mean = wsp_da.sel(time=slice(START_YEAR, END_YEAR)).mean("member").mean("time")
    wsp_mean.wsp.attrs["units"] = r"m s$^{-1}$"
    wsp_mean.wsp.attrs["long_name"] = "Mean wind speed"
    wsp_mean.wsp.attrs["description"] = "Wind speed at near surface (normally 10m)"
    wsp_mean.to_netcdf(os.path.join(_folder_name("mean"), "wsp.nc"))


def get_vars(var_list: List[str], regen_success_list=False):
    """
    Get the variable means and ensembles for the enviroment.

    Args:
        var_list (List[str]): List of variables to make.
    """
    for var_str in var_list:
        GetEnsemble(
            var=var_str,
            output_folder=_folder_name(var_str),
            regen_success_list=regen_success_list,
        )
        mean_var(var=var_str)


if __name__ == "__main__":
    # from src.data_loading.pangeo import GetEnsemble
    # python src/data_loading/pangeo.py
    # for var_str in VAR_PROP_D:
    #    GetEnsemble(
    #        var=var_str, output_folder=_folder_name(var_str), regen_success_list=True
    #    )
    # for var_str in VAR_PROP_D:
    #    mean_var(var=var_str)
    # make_wsp()
    # from src.data_loading.pangeo import get_vars
    get_vars(
        ["hurs", "psl"],
    )
