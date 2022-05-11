"""
ECMWF download script.
"""
import os
import xarray as xr
import cdsapi
from src.constants import DATA_PATH
from src.utils import timeit

DATA_DIREC = DATA_PATH / "ecmwf"


@timeit
def get_ecmwf(
    variable="total_precipitation",
    area=[
        90,
        -180,
        -90,
        180,
    ],
) -> None:
    """
    Get ECMWF variable.

    Args:
        variable (str, optional): ECMWF API variable name.
    Defaults to "total_precipitation".
    """

    if os.path.exists(DATA_DIREC):
        os.mkdir(DATA_DIREC)

    c = cdsapi.Client()

    c.retrieve(
        "reanalysis-era5-single-levels-monthly-means-preliminary-back-extension",
        {
            "format": "netcdf",
            "product_type": "reanalysis-monthly-means-of-daily-means",
            "variable": variable,
            "year": [
                "1950",
                "1951",
                "1952",
                "1953",
                "1954",
                "1955",
                "1956",
                "1957",
                "1958",
                "1959",
                "1960",
                "1961",
                "1962",
                "1963",
                "1964",
                "1965",
                "1966",
                "1967",
                "1968",
                "1969",
                "1970",
                "1971",
                "1972",
                "1973",
                "1974",
                "1975",
                "1976",
                "1977",
                "1978",
            ],
            "month": [
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
            ],
            "time": "00:00",
            "area": area,
        },
        str(DATA_DIREC / "1.nc"),
    )

    c.retrieve(
        "reanalysis-era5-single-levels-monthly-means",
        {
            "format": "netcdf",
            "product_type": "monthly_averaged_reanalysis",
            "variable": variable,
            "year": [
                "1979",
                "1980",
                "1981",
                "1982",
                "1983",
                "1984",
                "1985",
                "1986",
                "1987",
                "1988",
                "1989",
                "1990",
                "1991",
                "1992",
                "1993",
                "1994",
                "1995",
                "1996",
                "1997",
                "1998",
                "1999",
                "2000",
                "2001",
                "2002",
                "2003",
                "2004",
                "2005",
                "2006",
                "2007",
                "2008",
                "2009",
                "2010",
                "2011",
                "2012",
                "2013",
                "2014",
                "2015",
                "2016",
                "2017",
                "2018",
                "2019",
                "2020",
                "2021",
                "2022",
            ],
            "month": [
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
            ],
            "time": "00:00",
            "area": area,
        },
        str(DATA_DIREC / "2.nc"),
    )

    ds = xr.merge(
        [
            xr.open_dataset(str(DATA_DIREC / "1.nc")),
            xr.open_dataset(str(DATA_DIREC / "2.nc")),
        ]
    )
    ds.to_netcdf(DATA_DIREC / str(variable + ".nc"))
    os.remove(str(DATA_DIREC / "1.nc"))
    os.remove(str(DATA_DIREC / "2.nc"))


if __name__ == "__main__":
    # python src/data_loading/ecmwf.py
    get_ecmwf()
    main_variables = [
        "skin_temperature",
        "total_cloud_cover",
        "10m_wind_speed",
        "sea_surface_temperature",
        "mean_sea_level_pressure",
        "evaporation",
    ]

    other_variables = [
        "10m_u_component_of_wind",
        "10m_v_component_of_wind",
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
