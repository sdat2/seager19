"""File for formatting config.py file.


Example:
    Usage at top of program::

        from src.config.configs import format_config

    Instantiation somewhere in the progam::

        cfg = format_config(cfg)

"""
from omegaconf import DictConfig
from scipy.constants import pi, Stefan_Boltzmann, zero_Celsius, day
from astropy.constants import R_earth, g0


def derived_param(cfg: DictConfig) -> DictConfig:
    """
    Calculate dervied paramters using fundamental constants.

    Args:
        cfg (DictConfig): constants configuration.

    Returns:
        DictConfig: the config with all the derived param added.
    """
    cfg.atm["temp_surface_bar"] = zero_Celsius + cfg.atm.temp_surface_bar_celsius
    cfg.atm["qlh_coeff"] = cfg.atm.rho_air * cfg.atm.c_e * cfg.atm.latent_heat_vap
    cfg.atm["qlw_coeff"] = cfg.atm.emmisivity * Stefan_Boltzmann
    cfg.atm["eps_p"] = (
        pi / cfg.atm.height_tropopause ** 2 / cfg.atm.nbsq / cfg.atm.k_days / day
    )
    cfg.atm["eps"] = 1.0 / cfg.atm.eps_days / day
    cfg.atm["eps_u"] = cfg.atm.eps
    cfg.atm["eps_v"] = cfg.atm.eps * cfg.atm.e_frac
    cfg.atm["b_coeff"] = (
        g0.N * pi / cfg.atm["nbsq"] / cfg.atm["theta_00"] / cfg.atm["height_tropopause"]
    )
    cfg.atm["newtonian_cooling_coeff_k1"] = cfg.atm.b_coeff / cfg.atm.k_days / day
    cfg.atm["omega_2"] = 2 * (2 * pi / day)
    cfg.atm["pr_max"] = cfg.atm.pr_max_mm_day / day
    cfg.atm["beta"] = cfg.atm["omega_2"] / R_earth.N
    cfg.atm["y_south_lim"] = -cfg.atm.y_north_lim
    cfg.atm["dx"] = 360 / cfg.atm.nx
    cfg.atm["dy"] = (cfg.atm.y_north_lim - cfg.atm.y_south_lim) / cfg.atm.ny
    return cfg


def format_config(cfg: DictConfig) -> DictConfig:
    """cfg reformatting.

    Currently just evaluates the arithemtic put into the program.
    This isn't the best idea in the world, but it works, and
    it's broadly your own sorry fault if you put in
    dangerous code to the command line and it gets evaluated.

    Args:
        cfg (DictConfig): cfg to reformat.

    """

    # print("OmegaConf.to_yaml(cfg)", OmegaConf.to_yaml(cfg))
    # print(cfg.__repr__())

    return derived_param(cfg)
