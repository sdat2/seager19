"""
Zenodo download scripts.
"""
import requests
from tqdm import tqdm


# r = requests.get("https://sandbox.zenodo.org/record/1046378")
# print(r.status_code)
# from zenodo_get import zenodo_get
# import zenodopy
# https://sandbox.zenodo.org/record/1046378#.Ykx2f5NufmE
url = "https://sandbox.zenodo.org/record/1046378/files/CMIP6_INPUT_DATA.zip?download=1"

# zenodo_get(['param1','param2',....])
# zeno = zenodopy.Client()
# zeno.set_project('5867022')
# print(zeno.list_files)
# zeno.download_file(filename)
def get_zip() -> None:
    req = requests.get(url, stream=True)
    with open("CMIP6_INPUT_DATA.zip", "wb") as file:
        for chunk in tqdm(req.iter_content(chunk_size=128)):
            file.write(chunk)


if __name__ == "__main__":
    # python src/data_loading/zenodo.py
    get_zip()
    # print(r.status_code)
    # print(r.json())
    # print(zenodo_get.__doc__)
    # print(zeno.list_files)
    # print(zenodo_get(["-r=1046378", "-w" "-s"]))
