"""
ECMWF ERA5 download script.

This script requires you to have registered for an account at:

https://cds.climate.copernicus.eu/

And to have added the login details to your computer user profile as described there.

The script can take a long time to run, as requests at the CDS can
be held in a very long queue (i.e for many hours on weekdays).
"""
from typing import List, Tuple
import os
import shutil
import xarray as xr
import cdsapi
from src.constants import DATA_PATH, GWS_DIR
from src.utils import timeit

# from src.xr_utils import can_coords
# from src.data_loading.regrid import regrid_1d_to_standard


# intial directory on jasmin home workspace
DATA_DIREC = DATA_PATH / "ecmwf"
# final directory for archiving
ARCHIVE_DIREC = GWS_DIR / "ecmwf"

# regional bounding boxes
# Mekong river (some padding)
MEKONG = [35, 92, 8, 112]
# Gulf of Mexico
GOM = [-100, 15, 35, -80]


class FileNames:
    """
    Class to store ERA5 file names, both to create the files, and to
    reference them later.


    Properties:

    back_extension_path
    main_era5_path
    initial_combined_path
    archive_combined_path

    Example Usage::

        import xarray as xr
        from src.data_loading.ecmwf import FileNames
        precip_names = FileNames(variable="total_precipitation", region="mekong")
        xr.open_dataset(precip_names.archive_combined_path)


    """

    def __init__(self, variable: str = "total_precipitation", region: str = "") -> None:
        """
        Establish class with variable and region name.

        Args:
            variable (str, optional): Variable name.
        Defaults to "total_precipitation".
            region (str, optional): Region name. Defaults to "".
        """
        self.variable = variable
        # file paths
        if region != "":
            region = "_" + region
        self.back_extension_path = str(
            DATA_DIREC / str(variable + region + "_back_extension_era5.nc")
        )
        self.main_era5_path = str(DATA_DIREC / str(variable + region + "_main_era5.nc"))
        self.initial_combined_path = str(
            DATA_DIREC / str(variable + region + "_era5.nc")
        )
        self.archive_combined_path = str(
            ARCHIVE_DIREC / str(variable + region + "_era5.nc")
        )


def _year_lists(
    start_year: int = 1950, end_year: int = 2022, transition_year: int = 1979
) -> Tuple[List[str]]:
    """
    Return the year lists.

    Args:
        start_year (int, optional): First year in list. Defaults to 1950.
        end_year (int, optional): Final year in list. Defaults to 2022.
        transition_year (int, optional): Year to go between lists.
        Defaults to 1979. first year for newer ECMWF.


    Returns:
        Tuple[List[str]]: back_extension_years, main_era5_years

    Example::
        back_extension_years, main_era5_years = _year_lists(start_year=1951, end_year=2011)
    """
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
            for x in range(
                max(transition_year, start_year),
                min(end_year + 1, 2022 + 1)
                # 1
            )
        ]
    else:
        main_era5_years = []

    return back_extension_years, main_era5_years


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
    region: str = "",
    start_year: int = 1950,
    end_year: int = 2022,
    download: bool = True,  # whether to redownload data (or just archive)
    regrid: bool = False,  # whether to regrid onto 1 degree mesh.
    archive: bool = True,  # whether to archive in the group work space.
) -> None:
    """
    Get ECMWF monthly average variable.

    Args:
        variable (str, optional): ECMWF API variable name.
    Defaults to "total_precipitation".
        area (List[int], optional): Defaults to global [ 90, -180, -90, 180].
        region (str, optional): Region name to add to files. Defaults to "".
        start_year (int, optional): Start year of timeseries. Defaults to 1950.
        end_year (int, optional): End year of timeseries. Defaults to 2022.
        regrid (int, optional): Whether or not to regrid the data to my standard grid.
        archive (bool, optional): Defaults to True.
    """

    files = FileNames(variable=variable, region=region)

    ds_list = []

    back_extension_years, main_era5_years = _year_lists(
        start_year=start_year, end_year=end_year
    )

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
                files.back_extension_path,
            )
            ds_back = xr.open_dataset(files.back_extension_path)
            if "expver" in ds_back:
                ds_back.drop("expver")
            ds_list.append(ds_back)

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
                files.main_era5_path,
            )
            ds_main = xr.open_dataset(files.main_era5_path)
            print(ds_main)
            if "expver" in ds_main:
                ds_main = ds_main.isel(expver=0).drop("expver")
            ds_list.append(ds_main)
        # create new combined dataset.
        ds = xr.merge(ds_list)
        if regrid:
            print("Regridding not yet implemented")
            assert False
        ds.to_netcdf(files.initial_combined_path)
        # remove intermediate datasets.
        os.remove(files.back_extension_path)
        os.remove(files.main_era5_path)
        # move merged dataset to archive
    if archive:
        _archive_f(files)


def _archive_f(files: FileNames) -> None:
    """Archive file to gws.

    Args:
        files (FileNames): file structure class
    """
    shutil.move(
        files.initial_combined_path,
        files.archive_combined_path,
    )


def _dearchive_f(files: FileNames) -> None:
    """Dearchive file from gws.

    Args:
        files (FileNames): file structure class.
    """
    shutil.move(
        files.archive_combined_path,
        files.initial_combined_path,
    )


@timeit
def get_main_variables(regrid=False) -> None:
    """Make the main ERA5 variables for seager19."""
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
        get_era5(variable=var, regrid=regrid)


MEKONG_VARIABLES = [
    "evaporation",
    "total_precipitation",
    "skin_temperature",
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


@timeit
def get_mekong_variables() -> None:
    """Download and archive mekong variables"""
    for var in MEKONG_VARIABLES:
        get_era5(variable=var, area=MEKONG, region="mekong")


@timeit
def rename_mekong() -> None:
    """Rename mekong files."""
    for var in MEKONG_VARIABLES:
        wrong_names = FileNames(variable=var, region="mekong_")
        right_names = FileNames(variable=var, region="mekong")
        shutil.move(
            wrong_names.archive_combined_path, right_names.archive_combined_path
        )


def _test_year_lists() -> None:
    """Test if year lists output correct results."""
    assert str(_year_lists(1978, 1978)) == "(['1978'], [])"
    assert str(_year_lists(1979, 1979)) == "([], ['1979'])"


if __name__ == "__main__":
    # python src/data_loading/ecmwf.py
    rename_mekong()
    # get_main_variables()
    # get_mekong_variables()
    # _test_year_lists()
    # print(_year_lists(1980, 2011))
    # print(_year_lists(1960, 1970))
    # print(_year_lists(1950, 2022))
    # print(_year_lists(1978, 1978))
    # print(_year_lists(1979, 1979))
    # print(
    #    FileNames(
    #        variable="total_precipitation", region="mekong"
    #    ).archive_combined_path
    # )

    # get_era5(variable="skin_temperature", download=True)
