"""Make fortran90 documentation?"""
from sphinxfortran import fortran_autodoc
import os


# list_of_src_files = [
# "ocean/SRC/data-mod.F",
# "ocean/SRC/om_ekm.F",
# "ocean/SRC/om_equi.F",
# "om_equi.F",
# "ocean/SRC/om_forc.F",
# "ocean/SRC/om_leap.F",
# "ocean/SRC/om_main.F",
# "ocean/SRC/om_mem.F",
# "ocean/SRC/om_sst.F",
# "ocean/SRC/om_qflux.F",
# "ocean/SRC/om_tios.F",
# "ocean/SRC/om_wrap.F",
# "ocean/SRC/wrap-mod.F",
# ]

"""
list_of_src_files = [
    "data-mod.F",
    "om_ekm.F",
    "om_equi.F",
    "om_forc.F",
    "om_leap.F",
    "om_main.F",
    "om_mem.F",
    "om_sst.F",
    "om_qflux.F",
    "om_tios.F",
    "om_wrap.F",
    "wrap-mod.F",
]
"""

list_of_src_files = [x for x in os.listdir() if x.endswith(".F")]

# list_of_src_files = ["../" + x for x in list_of_src_files]

print(list_of_src_files)

fortran_autodoc.F90toRst(list_of_src_files)
