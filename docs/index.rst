.. src documentation master file, created by
   sphinx-quickstart on Thu Mar 18 17:37:46 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Analysis of a Parsimonious Coupled Model of the Equatorial Pacific Surface Temperature Change
==============================

Welcome to the seager19 documentation!

This project reassembles a parsimonious coupled model of the equatorial Pacific, 
from Seager et al. 2019 (S19), that was created to explain the cold 
tongue bias in CMIP5 models. When forced with ECMWF reanalysis fields, it can reproduce
the trend observed in ECMWF/ORAS4 reanalysis product that was forced with the same
fields. It shows that the CMIP5 bias inthe trend in NINO3.4 from 1958-2017 could
be due to a product of the CMIP5 bias in relative humidityand sea surface winds, 
which is shown through exchanging ECMWF mean fields for CMIP5 multimodelmean fields.
The replacements of mean relative humidity, mean wind speed, and both together, 
lead toincreases in the NINO3.4 trend of 0.31±0.03 K, 0.054±0.005 K, and 0.47±0.04 K 
respectively when testedwith a range of plausible inputs. This is congruent with the
observed difference of 0.478 K between theECMWF/ORAS4 reanalysis product and the CMIP5 multimodel mean. I investigate how reliable the results from this model might be by varying
the free parameters and findthat, as far as tested, the model is not overly sensitive to
subjective inputs. It is therefore plausible thatobserved bias in the increase in sea
surface temperature is caused by excess humidity, and insufficient windspeed over the cold
tongue – are reinforced as credible, as the model is not overly sensitive to variation 
of the free parameters.

The first section `seager19` contains the main `README.md` of the repository,
so as to reduce duplication, and this should provide a reasonable introduction
to the repository as a whole.

Here is the current breakdown of the model code by language:

.. code-block:: bash

   $ cloc --report-file=docs/lang.txt  $(git ls-files)

.. include:: lang.txt
   :literal:


MRes Proposal (See :download:`the final report <Report_without_documentation.pdf>` for results.):

Seager et al. 2019 [1, hereafter S19] showed that although CMIP5 ensemble members have
a positive NINO3.4 trend (towards El Nino) where as the observations show a more 
negative NINO3.4 trend (towards La Nina).

.. image:: gifs/trend_graph.png
  :width: 500
  :alt: trend graph

Caption: observations = orange diamonds / blue crosses;
models = blue / black. As you can see the CMIP ensemble 
members show little overlap with any of the renalysis products.

.. image:: gifs/trend_graphic.png
  :width: 500
  :alt: trend graphic

Caption: This suggests an over all tendency to La Nina in observations rather than El Nino.


They showed that the observed trend can be reproduced with a simple coupled physical model.
Here, we carry out a parameter sensitivity analysis of the S19 model. 
Of particular interest might be the S19 model's sensitivity to the drag coefficient, 
as S19 note that they chose a much higher value than normal so as to replicate the
amplitude of ENSO. This sensitivity analysis could first be achieved using a Gaussian
Process (GP) with a radial basis function (RBF) kernel of a given smoothness,
as the number of data points will initially be quite small (<10^{4}).
S19 is computationally lightweight, allowing for a large number of parallel
sensitivity experiments to be run at the same time in order to generate the
training dataset for the GP model. The GP model will allow us to rapidly 
explore the parameter space in between our chosen parameter configurations, 
in terms of both the mean value and uncertainties. From this initial baseline, 
we could expand to more sophisticated sensitivity analyses, and/or more complicated 
model settings.


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
