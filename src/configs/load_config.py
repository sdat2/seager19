"""Load hydra config (util for tests)."""
from hydra.experimental import initialize, compose
from src.constants import PROJECT_PATH, CONFIG_NAME, CONFIG_PATH
from src.configs.config import format_config
from omegaconf import DictConfig


# 1. initialize will add config_path the config search path within the context
# 2. The module with your configs should be importable.
#    it needs to have a __init__.py (can be empty). pytest src/test/test_hydra.py
# 3. THe config path is relative to the file calling initialize (this file)
def load_config(prefix: str = "../../") -> DictConfig:
    """Tests loading hydra config file."""
    with initialize(
        config_path=prefix + str(CONFIG_PATH).replace(str(PROJECT_PATH) + "/", "")
    ):
        # config is relative to a module
        cfg = compose(
            config_name=CONFIG_NAME, overrides=["user=test_user", "name=test_run"]
        )
        print(cfg)
        cfg = format_config(cfg)
        print(cfg)

    return cfg
