.. src documentation master file, created by
   sphinx-quickstart on Thu Mar 18 17:37:46 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Home - seager19 Documentation
==============================
Welcome to the seager19 documentation!

This webpage/document explores the code and data 
used in the seager19 MRes project, and some of the 
motivations for different choices made. The 
first section `seager19` contains the main
`README.md` of the repository, so as to reduce duplication,
and this should provide a reasonable introduction to the 
repository as a whole.

The online version of this documentation is able to include
more interactive features such as gifs of the results.
Therefore, I suggest that the web version is 
probably a better read.


MRes Proposal
-------------

Seager et al. 2019 [1, hereafter S19] showed that although CMIP5 ensemble members have a positive NINO3.4 trend where as the observations show a negative NINO3.4 trend.

.. image:: gifs/trend_graph.png
  :width: 500
  :alt: trend graph

Caption: observations = orange diamonds / blue crosses;
 models = blue / black.

.. image:: gifs/trend_graphic.png
  :width: 500
  :alt: trend graphic

Caption: this suggests an over all tendency to La Nina in observations rather than El Nino.

They showed that the observed trend can be reproduced with a simple coupled physical model.
Here, we carry out a parameter sensitivity analysis of the S19 model. Of particular interest might be the S19 model's sensitivity to the drag coefficient, as S19 note that they chose a much higher value than normal so as to replicate the amplitude of ENSO. This sensitivity analysis could first be achieved using a Gaussian Process (GP) with a radial basis function (RBF) kernel of a given smoothness, as the number of data points will initially be quite small (<10^{4}). S19 is computationally lightweight, allowing for a large number of parallel sensitivity experiments to be run at the same time in order to generate the training dataset for the GP model. The GP model will allow us to rapidly explore the parameter space in between our chosen parameter configurations, in terms of both the mean value and uncertainties. From this initial baseline, we could expand to more sophisticated sensitivity analyses, and/or more complicated model settings.


.. image:: gifs/om_diag_SST_SST.gif
  :width: 500
  :alt: SST in setup period

Scientific questions to be addressed include:

- Can we replicate the results displayed in S19?
- How robust is the model to the parameters chosen?
- Can the sensitivity of the model to the parameters be understood from the physical processes underlying it?
- How skillful are different emulation functions at fitting input/output of the model?

Skills learnt
-------------

This project will be a good opportunity to learn about numerical solutions to differential equations [3, 4]. From the rubric at https://ai4er-cdt.esc.cam.ac.uk this project would include the primary application areas of Weather, Climate and Air Quality, as well as touching on _Natural Hazards_ as this bias in CMIP will effect the frequency and intensity of tropical cyclones (and therefore their effects such as storm surges, as in my MSci thesis https://bit.ly/msci-report).  The methodology of this project will primarily be within 2. Environmental modelling.  The model used is quite low down in the hierarchy of models, and so it will be much easier to demonstrate simple ML tools upon it.

MRes Citations
--------------

[1] Seager,  R. et  al.  Strengthening  tropical  Pacific  zonal  sea  surface  temperature  gradient  consistent  with rising  greenhouse  gases, July 2019. https://doi.org/10.1038/s41558-019-0505-x

[2] Tian, B. & Dong, X. The Double-ITCZ Bias in CMIP3, CMIP5, and CMIP6 Models Based on Annual Mean Precipitation. Geophysical  Research  Letters 47, e2020GL087232. issn: 0094-8276. doi:10.1029/2020GL087232. https://onlinelibrary.wiley.com/doi/abs/10.1029/2020GL087232 (Apr. 2020).

[3] Hinch,  E. J. Think  Before  You  Compute doi:10.1017/9781108855297. https://www.cambridge.org/core/books/think-before-you-compute/8CC661846951CE3F08F1CC739FBB341F (Cambridge University Press, Apr. 2020).

[4] Iserles, A. Numerical  Solution  of  Differential  Equations Lent Term 2021 (Centre for Mathematical Sciences, University of Cambridge, 2021)


.. toctree::
   :maxdepth: 3
   :caption: Contents:

   self
   MAIN_README.md
   paper-fig.ipynb
   src
   current
   gallery.md
   dQ.ipynb
   TCAM.ipynb
   about


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
