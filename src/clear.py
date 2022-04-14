"""Script for clearing all the finished runs from the SSD."""
import os
from src.wandb_utils import finished_names
from src.constants import LOG_PATH, DEFAULT_PROJECT


def clear(project: str = DEFAULT_PROJECT) -> None:
    """
    Clear the logs.

    Args:
        project (str): Wandb project. Defaults to src.constants.DEFAULT_PROJECT.
    """
    rem_list = [i for i in finished_names(project=project) if i in os.listdir(LOG_PATH)]
    print(rem_list)
    for i in rem_list:
        # delete_folder(os.path.join(LOG_PATH, i))
        os.system("rm -rf " + os.path.join(LOG_PATH, i))


if __name__ == "__main__":
    # python src/clear.py
    clear(project="sdat2/ENSOTrend-beta")
