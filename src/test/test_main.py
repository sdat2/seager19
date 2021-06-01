"""test main."""
from src.main import sub_main
from src.configs.load_config import load_config


def test_main() -> None:
    """Test main."""
    cfg = load_config()
    sub_main(cfg, unit_test=True)
