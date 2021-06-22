"""Test wandb."""
from src.wandb_utils import start_wandb
from src.configs.load_config import load_config


def test_wandb_util() -> None:
    """Test src.configs.wandb_util"""
    cfg = load_config()
    start_wandb(cfg, unit_test=True)
