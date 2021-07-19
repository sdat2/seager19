"""Script for clearing all the finished runs from the SSD."""
import os
from src.wandb_utils import finished_names
from src.constants import LOG_PATH

if __name__ == "__main__":
    # python src/clear.py
    rem_list = [i for i in finished_names() if i in os.listdir(LOG_PATH)]
    print(rem_list)
    for i in rem_list:
        # delete_folder(os.path.join(LOG_PATH, i))
        os.system("rm -rf " + os.path.join(LOG_PATH, i))
