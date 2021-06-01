"""Test ocean model runs."""
from src.configs.load_config import load_config
from src.models.model_setup import ModelSetup
from src.constants import TEST_DIREC
from src.models.ocean import Ocean


def test_ocean() -> None:
    """Test ocean model runs."""
    cfg = load_config()
    setup = ModelSetup(str(TEST_DIREC))
    ocean = Ocean(cfg, setup)
    ocean.compile_all()
    if cfg.run:
        ocean.run_all()
    if cfg.animate:
        ocean.animate_all()
