"""Pangeo.py

Possible preprocessing script for pangeo config script.

# TODO: create auto-archive functionality.

TODO: change Filestructure to be more logical.

File structures inside src/data:

Initial directory:
    DIR=src/data/nc

Final directory:
    DIR=/gws/nopw/j04/ai4er/users/sdat2/CMIP6

    Full ensemble time series:
        ${DIR}/historical.ssp585/timeseries/${var}/${member}.nc

    Ensemble 60 year mean state:
        ${DIR}/historical.ssp585/mean/${var}/${member}.nc

    Ensemble 60 year climatology:
        ${DIR}/historical.ssp585/climatology/${var}/${member}.nc

    Ensemble 60 year linear trend:
        ${DIR}/historical.ssp585/trend/${var}/${member}.nc

    Ensemble 60 year linear trend:
        ${DIR}/historical.ssp585/nino3.4/${var}/${member}.nc

    Ensemble 60 year nino3.4:
        ${DIR}/historical.ssp585.nino3.4/${var}.nc

    Ensemble 60 year nino3.4.climatology:
        ${DIR}/historical.ssp585.nino3.4.climatology/${var}.nc

    Ensemble 60 year nino3.4.trends:
        ${DIR}/historical.ssp585.nino3.4.trends/${var}.nc

    MMM 60 year mean:
        ${DIR}/historical.ssp585.mmm/mean/${var}.nc

    MMM 60 year trend:
        ${DIR}/historical.ssp585.mmm/trend/${var}.nc

    MMM 60 year climatology:
        ${DIR}/historical.ssp585.mmm/climatology/${var}.nc

"""
