"""Load hydra config (util for tests)."""
from hydra import initialize, compose  # new hydra moved from experimental
from src.constants import PROJECT_PATH, CONFIG_NAME, CONFIG_PATH
from src.configs.config import format_config
from omegaconf import DictConfig


# 1. initialize will add config_path the config search path within the context
# 2. The module with your configs should be importable.
#    it needs to have a __init__.py (can be empty). pytest src/test/test_hydra.py
# 3. THe config path is relative to the file calling initialize (this file)
def load_config(prefix: str = "../../", test=True) -> DictConfig:
    """
    Tests loading hydra config file.

    Args:
        prefix (str, optional): . Defaults to "../../".
        test (bool, optional): Whether on not this is a test.
            Defaults to True. This will turn off running the
            main sections of the ocean model.

    Returns:
        DictConfig: Return the dict config.
    """
    with initialize(
        config_path=prefix + str(CONFIG_PATH).replace(str(PROJECT_PATH) + "/", "")
    ):
        # config is relative to a module
        if test:
            override_list = [
                "user=test_user",
                "name=test_run",
                "ocean.spin=false",
                "ocean.diag=false",
                "ocean.ingrid=false",
                "ocean.run_through=false",
            ]
        else:
            override_list = ["user=non-test"]
        cfg = compose(
            config_name=CONFIG_NAME,
            overrides=override_list,
        )
        # print(cfg)
        cfg = format_config(cfg)
        # print(cfg)

    return cfg
