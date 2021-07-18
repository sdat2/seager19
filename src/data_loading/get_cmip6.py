"""Get CMIP6 variables on pangeo."""
import os
from typing import Union
import numpy as np
import pandas as pd
import xarray as xr
import dask
from cmip6_preprocessing.preprocessing import (
    combined_preprocessing,
)
import cftime
from intake import open_catalog


@np.vectorize
def standardise_time(
    time: Union[cftime.datetime, np.datetime64],
    calendar="standard",  # "gregorian"
) -> cftime.datetime:
    """
    Standardise time.

    Args:
        time (Union[ cftime._cftime.DatetimeNoLeap, cftime._cftime.Datetime360Day,
                    np.datetime64 ]): Time array.
        calendar (str, optional): Which cftime calendar to replace it with.
            Defaults to "standard".

    Returns:
        cftime._cftime.Datetime360Day: The new calendar.
    """
    if isinstance(time, np.datetime64):
        time = pd.to_datetime(time)
    # put the new time in the middle of the given month
    # return np.datetime64(str(str(time.year) +
    #  '-' + str(time.month) + '-' + str(15) + "T00:00:00"))
    return cftime.datetime(time.year, time.month, 15, calendar=calendar)  # "360_day")


def preproc(ds: Union[xr.Dataset, xr.DataArray]) -> Union[xr.Dataset, xr.DataArray]:
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
        if not os.path.exists(output_folder):
            os.mkdir(output_folder)

        self.var = var
        self.cat = open_catalog(
            str(
                "https://raw.githubusercontent.com/pangeo-data/"
                + "pangeo-datastore/master/intake-catalogs/master.yaml"
            )
        )["climate"]["cmip6_gcs"]
        self.instit = self.cat.unique(["institution_id"])["institution_id"]["values"]
        default_success_list = [
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
        if regen_success_list:
            self.success_list = self.get_sucess_list()
        else:
            self.success_list = default_success_list
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
                        zarr_kwargs=z_kwargs, preprocess=preproc
                    )

                # print(i, "suceeded")
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
        experiment="historical",
        year_begin="1958",
        year_end="2014",
        var="ts",
        # pylint: disable=dangerous-default-value
        xlim=[100, 290],
        ylim=[-30, 30],
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
                x
                for x in self.success_list
                if x not in ["AWI", "MRI", "CSIRO-ARCCSS", "CCCma"]
            ],
        )
        subset = self.cat.search(**query)

        z_kwargs = {"consolidated": True, "decode_times": True}

        # pass the preprocessing directly
        with dask.config.set(**{"array.slicing.split_large_chunks": True}):
            dset_dict_proc = subset.to_dataset_dict(
                zarr_kwargs=z_kwargs, preprocess=preproc
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
                if "height" in [c for c in sub_da.coords]:
                    sub_da = sub_da.drop("height")
                da_list.append(sub_da)
        if self.var in ["uas", "vas"]:
            da = xr.concat(da_list, "member", coords="minimal")
        else:
            da = xr.concat(da_list, "member")
        da = da.assign_coords({"member": key_list})  # , "time": times})

        return da

    def get_ensemble(self, var: str = "ts") -> xr.DataArray:
        """
        Get ts.

        Args:
            var (str, optional): The variable. Defaults to "ts".

        Returns:
            xr.DataArray: Ensemble datarray full of different items.
        """
        da_hist = self.get_var(xlim=[0, 360], ylim=[-80, 80], var=var)
        da_ssp585 = self.get_var(
            experiment="ssp585",
            year_begin="2014",
            year_end="2017",
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
                            + str(instit)
                            + "."
                            + str(model)
                            + "."
                            + str(member_id)
                            + ".80.nc",
                        )
                    )
                else:
                    print("problem with " + model + " " + member_id)
        # da = xr.concat(da_list, "member")
        # da = da.assign_coords({"member": key_list})
        # return da_list


if __name__ == "__main__":
    # from src/data_loading/get_cmip6 import GetEnsemble
    GetEnsemble(var="ps", output_folder="nc_ps")
