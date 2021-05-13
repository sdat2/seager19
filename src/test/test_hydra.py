"""Test hydra."""
from hydra.experimental import initialize, compose
from src.constants import CONFIG_NAME, CONFIG_PATH
from src.configs.config import format_config

# 1. initialize will add config_path the config search path within the context
# 2. The module with your configs should be importable.
#    it needs to have a __init__.py (can be empty).
# 3. THe config path is relative to the file calling initialize (this file)
def test_with_initialize() -> None:
    """Tests loading hydra config file."""
    with initialize(config_path=CONFIG_PATH):
        # config is relative to a module
        cfg = compose(config_name=CONFIG_NAME, overrides=["user=test_user"])
        print(cfg)
        cfg = format_config(cfg)
        print(cfg)
