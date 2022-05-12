"""
ECMWF ERA5 download script.
"""
from typing import List
import os
import shutil
import xarray as xr
import cdsapi
from src.constants import DATA_PATH, GWS_DIR
from src.utils import timeit

# intial directory on jasmin home workspace
DATA_DIREC = DATA_PATH / "ecmwf"
# final directory for archiving
ARCHIVE_DIREC = GWS_DIR / "ecmwf"


@timeit
def get_era5(
    variable: str = "total_precipitation",
    # pylint: disable=dangerous-default-value
    area: List[int] = [
        90,
        -180,
        -90,
        180,
    ],
    start_year: int = 1950,
    end_year: int = 2022,
    download: bool = True,  # whether to redownload data (or just archive)
    regrid: bool = False,  # whether to regrid onto 1 degree mesh.
    archive: bool = True,
) -> None:
    """
    Get ECMWF variable.

    Args:
        variable (str, optional): ECMWF API variable name.
    Defaults to "total_precipitation".
        area (List[int], optional): Defaults to [ 90, -180, -90, 180, ].
        archive (bool, optional): Defaults to True.
    """
    transition_year = 1979  # first year for newer ECMWF
    # make year lists splitting years around 1950
    if start_year < transition_year:
        back_extension_years = [
            str(x)
            for x in range(max(1950, start_year), min(end_year + 1, transition_year))
        ]
    else:
        back_extension_years = []
    if end_year >= transition_year:
        main_era5_years = [
            str(x)
            for x in range(max(transition_year, start_year), min(end_year + 1, 2022+1))
        ]
    else:
        main_era5_years = []

    back_extension_path = str(DATA_DIREC / str(variable + "1.nc"))
    main_era5_path = str(DATA_DIREC / str(variable + "1.nc"))
    initial_combined_path = str(DATA_DIREC / str(variable + ".nc"))
    archive_combined_path = str(ARCHIVE_DIREC / str(variable + ".nc"))

    ds_list = []

    print(back_extension_years, main_era5_years)

    # make sure the right directories exist.
    if not os.path.exists(DATA_DIREC):
        os.mkdir(DATA_DIREC)
    if archive:
        if not os.path.exists(ARCHIVE_DIREC):
            os.mkdir(ARCHIVE_DIREC)

    months = [
        "01",
        "02",
        "03",
        "04",
        "05",
        "06",
        "07",
        "08",
        "09",
        "10",
        "11",
        "12",
    ]

    if download:

        # Download data
        c = cdsapi.Client()

        if back_extension_years:  # wont run if empty
            c.retrieve(
                "reanalysis-era5-single-levels-monthly-means-preliminary-back-extension",
                {
                    "format": "netcdf",
                    "product_type": "reanalysis-monthly-means-of-daily-means",
                    "variable": variable,
                    "year": back_extension_years,
                    "month": months,
                    "time": "00:00",
                    "area": area,
                },
                back_extension_path,
            )
            ds_list.append(xr.open_dataset(back_extension_path).drop("expver"))
        if main_era5_years:  # wont run if empty
            c.retrieve(
                "reanalysis-era5-single-levels-monthly-means",
                {
                    "format": "netcdf",
                    "product_type": "monthly_averaged_reanalysis",
                    "variable": variable,
                    "year": main_era5_years,
                    "month": months,
                    "time": "00:00",
                    "area": area,
                },
                main_era5_path,
            )
            ds_list.append(xr.open_dataset(main_era5_path).drop("expver"))
        # create new combined dataset.
        ds = xr.merge(ds_list)
        if regrid:
            print("Regridding not yet implemented")
            assert False
        ds.to_netcdf(initial_combined_path)
        # remove intermediate datasets.
        os.remove(back_extension_path)
        os.remove(main_era5_path)
        # move merged dataset to archive
    if archive:
        shutil.move(
            initial_combined_path,
            archive_combined_path,
        )


@timeit
def get_main_variables():
    main_variables = [
        "total_precipitation",
        "skin_temperature",
        "total_cloud_cover",
        "10m_wind_speed",
        "sea_surface_temperature",
        "mean_sea_level_pressure",
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
    ]
    for var in main_variables:
        get_era5(variable=var)


if __name__ == "__main__":
    # python src/data_loading/ecmwf.py

    get_era5(variable="skin_temperature", download=True)

    other_variables = [
        "evaporation",
        "low_vegetation_cover",
        "soil_temperature_level_1",
        "soil_temperature_level_2",
        "soil_temperature_level_3",
        "soil_temperature_level_4",
        "soil_type",
        "type_of_high_vegetation",
        "volumetric_soil_water_layer_1",
        "volumetric_soil_water_layer_2",
        "volumetric_soil_water_layer_3",
        "volumetric_soil_water_layer_4",
    ]
