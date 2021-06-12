"""Move old logs."""
# import os
import shutil


def delete_folder(folder: str = "/path/to/folder") -> None:
    """
    Delete the folder.

    Args:
        folder (str, optional): Path to a folder. Defaults to "/path/to/folder".
    """
    try:
        shutil.rmtree(folder)
    except OSError as e:
        print("Error: %s - %s." % (e.filename, e.strerror))
