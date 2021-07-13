.. src documentation master file, created by
   sphinx-quickstart on Thu Mar 18 17:37:46 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Home - seager19 Documentation
==============================

Welcome to the seager19 documentation!

This webpage/document explores the code and data 
from Seager et al. 2019 (S19), reeassembling their coupled 
model of the tropical pacific. The 
first section `seager19` contains the main
`README.md` of the repository, so as to reduce duplication,
and this should provide a reasonable introduction to the 
repository as a whole.

The online version of this documentation is able to include
more interactive features such as gifs of the results.
Therefore, I suggest that the web version is 
probably a better read.

Current state of code:

.. code-block:: bash

   $ cloc --report-file=docs/lang.txt --sum-one  $(git ls-files)


.. include:: lang.txt
   :literal:


MRes Proposal (See :download:`this example script <Report_without_documentation>`.):

Seager et al. 2019 [1, hereafter S19] showed that although CMIP5 ensemble members have a positive NINO3.4 trend (towards El Nino) where as the observations show a more negative NINO3.4 trend (towards La Nina).

.. image:: gifs/trend_graph.png
  :width: 500
  :alt: trend graph

Caption: observations = orange diamonds / blue crosses;
 models = blue / black. As you can see the CMIP ensemble members show little overlap with any of the renalysis products.

.. image:: gifs/trend_graphic.png
  :width: 500
  :alt: trend graphic

Caption: this suggests an over all tendency to La Nina in observations rather than El Nino.

They showed that the observed trend can be reproduced with a simple coupled physical model.
Here, we carry out a parameter sensitivity analysis of the S19 model. Of particular interest might be the S19 model's sensitivity to the drag coefficient, as S19 note that they chose a much higher value than normal so as to replicate the amplitude of ENSO. This sensitivity analysis could first be achieved using a Gaussian Process (GP) with a radial basis function (RBF) kernel of a given smoothness, as the number of data points will initially be quite small (<10^{4}). S19 is computationally lightweight, allowing for a large number of parallel sensitivity experiments to be run at the same time in order to generate the training dataset for the GP model. The GP model will allow us to rapidly explore the parameter space in between our chosen parameter configurations, in terms of both the mean value and uncertainties. From this initial baseline, we could expand to more sophisticated sensitivity analyses, and/or more complicated model settings.


.. image:: gifs/om_diag_SST_SST.gif
  :width: 500
  :alt: SST in diag period 1956-58

Scientific questions to be addressed include:

- Can we replicate the results displayed in S19? [partially]
- How robust is the model to the parameters chosen? [fairly]
- Can the sensitivity of the model to the parameters be understood from the physical processes underlying it? [partially]
- How skillful are different emulation functions at fitting input/output of the model? [untested]

Citations:

[1] Seager,  R. et  al.  Strengthening  tropical  Pacific  zonal  sea  surface  temperature  gradient  consistent  with rising  greenhouse  gases, July 2019. https://doi.org/10.1038/s41558-019-0505-x

[2] Tian, B. & Dong, X. The Double-ITCZ Bias in CMIP3, CMIP5, and CMIP6 Models Based on Annual Mean Precipitation. Geophysical  Research  Letters 47, e2020GL087232. issn: 0094-8276. doi:10.1029/2020GL087232. https://onlinelibrary.wiley.com/doi/abs/10.1029/2020GL087232 (Apr. 2020).



.. toctree::
   :maxdepth: 3
   :caption: Contents:

   self
   MAIN_README.md
   src
   OCEAN_README.md
   gallery.md
   seager-rep.ipynb
   about


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
