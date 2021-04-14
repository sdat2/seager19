# MRes Proposal

![SST output over spin up period](gifs/om_diag_SST_SST.gif)

Seager et al. 2019 [1, hereafter S19] showed that although CMIP5 ensemble members do not follow the observed trend in ENSO, this observed trend can be reproduced with a simple coupled physical model. Here, we carry out a parameter sensitivity analysis of the S19 model. Of particular interest might be the S19 model's sensitivity to the drag coefficient, as S19 note that they chose a much higher value than normal so as to replicate the amplitude of ENSO. This sensitivity analysis could first be achieved using a Gaussian Process (GP) with a radial basis function (RBF) kernel of a given smoothness, as the number of data points will initially be quite small (<10^{4}). S19 is computationally lightweight, allowing for a large number of parallel sensitivity experiments to be run at the same time in order to generate the training dataset for the GP model. The GP model will allow us to rapidly explore the parameter space in between our chosen parameter configurations, in terms of both the mean value and uncertainties. From this initial baseline, we could expand to more sophisticated sensitivity analyses, and/or more complicated model settings.

Scientific questions to be addressed include:

- Can we replicate the results displayed in S19?
- How robust is the model to the parameters chosen?
- Can the sensitivity of the model to the parameters be understood from the physical processes underlying it?
- How skillful are different emulation functions at fitting input/output of the model?

## Skills learnt

This project will be a good opportunity to learn about numerical solutions to differential equations [3, 4]. From the rubric at <https://ai4er-cdt.esc.cam.ac.uk> this project would include the primary application areas of _Weather, Climate and Air Quality_, as well as touching on _Natural Hazards_ as this bias in CMIP will effect the frequency and intensity of tropical cyclones (and therefore their effects such as storm surges, as in my MSci thesis <https://bit.ly/msci-report>).  The methodology of this project will primarily be within _2. Environmental modelling_.  The model used is quite low down in the hierarchy of models, and so it will be much easier to demonstrate simple ML tools upon it.

### MRes Citations

[1] Seager,  R. et  al.  Strengthening  tropical  Pacific  zonal  sea  surface  temperature  gradient  consistent  with rising  greenhouse  gases, July 2019. <https://doi.org/10.1038/s41558-019-0505-x>.

[2] Tian, B. & Dong, X. The Double-ITCZ Bias in CMIP3, CMIP5, and CMIP6 Models Based on Annual Mean Precipitation. Geophysical  Research  Letters 47, e2020GL087232. issn: 0094-8276. doi:10.1029/2020GL087232. <https://onlinelibrary.wiley.com/doi/abs/10.1029/2020GL087232> (Apr. 2020).

[3] Hinch,  E. J. Think  Before  You  Compute doi:10.1017/9781108855297. <https://www.cambridge.org/core/books/think-before-you-compute/8CC661846951CE3F08F1CC739FBB341F> (Cambridge University Press, Apr. 2020).

[4] Iserles, A. Numerical  Solution  of  Differential  Equations Lent Term 2021 (Centre for Mathematical Sciences, University of Cambridge, 2021)

![Hurricane Katrina at Landfall: Credit NASA](https://cdn.britannica.com/74/121674-050-C458B2B5/satellite-image-National-Oceanic-and-Atmospheric-Administration-August-28-2005.jpg)

## Potential PhD project

### Hybrid tropical cyclones hazard modelling, focussing on storm surges on the US East coast

Tropical cyclones (TCs) have historically been the world’s largest physical hazard [2] in terms of economic damage and specifically TC storm surges for lives lost [1, 2, 3]. This hazard is increasing as more people live on vulnerable coastlines in substandard buildings [1]. Anthropogenic climate change is expected to exacerbate this, as it begins to increase the maximum TC intensity and range [4, 5, 6, 7].

It is an open question whether there are ways of summarising the
characteristics of a tropical cyclone and the coastline so that good statistical models can be built to assess their impact, as a more efficient alternative to high resolution 3D modelling. Particular features, such as the distance to the 30m contour have been used as informative statistical features for hazard models (see e.g. [10]).
My MSci Project <https://bit.ly/msci-report>, <https://bit.ly/msci-summary>
showed that you can predict the responsiveness of a coastline to a windstess fairly well with a simple physical model:

![Responsiveness metric: for further details see msci summary.](gifs/responsiveness.png)

This metric captured the broad trends seen in the model output. As a possible improvement, we propose to develop better feature engineering to incorporate elements such as the convexity of the coast, in a model similar to Chavas et al. [10]. Talea Mayo (Emory University) has offered to provide an initial storm surge dataset that was produced in collaboration with Ning Lin [11, 12, 13], and this could provide an initial dataset to address the scientific question:

- Are there metrics that can summarise the bathymetry’s impact on the size of the storm surge from a given tropical cyclone?
- More generally, is it possible to build a statistical or machine learning type model to connect a variety of features (e.g. bathymetry, local winds) to storm surge risk?
- Even more generally, can we extend this framework to other forms of risk (e.g. wind damage)?

## Additional analysis of CMIP bias

A specific objective in the early part of the project would be to use a cloud computing platform (e.g. JASMIN’s implementation of PANGEO) to examine CMIP6 for the ENSO bias discussed in Seager et al. (2019). Initial results indicate that this bias does persist in CMIP6 (e.g. [11]).

Scientific questions include:

- Do any particular model families perform better, and is there a physical reason for this?
- What are the likely effects of these biases on risk from e.g. Tropical cyclones (using a metric like the potential intensity)?

### PhD Citations

[1] Emanuel, K. et al. Divine Wind: The History  and Science of Hurricanes (Oxford university press, 2005).

[2] Shultz, J. M., Russell, J. & Espinel, Z. Epidemiology of tropical cyclones: the dynamics of disaster, disease, and development. Epidemiologic reviews 27, 21–35 (2005).

[3] Zhang, Q., Wu, L. & Liu, Q. Tropical cyclone damages in China 1983–2006. Bulletin of the American Meteorological Society 90, 489–496 (2009).

[4] Mendelsohn, R., Emanuel, K. A., Chonabayashi, S. & Bakkensen, L. The impact of climate change on global tropical cyclone damage. Nature climate change 2, 205–209 (2012).

[5] Emanuel, K. A., Sundararajan, R. & Williams, J. Hurricanes and global warming: Results from downscaling IPCC AR4 simulations. Bulletin of the American Meteorological Society 89, 347–368 (2008).

[6] Emanuel, K. A. Will global warming make hurricane forecasting more di
cult? Bulletin of the American Meteorological Society 98, 495–501 (2017).

[7] NORDHAUS, W. D. The economics of hurricanes and implications of global warming. Climate Change Economics 01, 1–20 (2010).

[8] Fedorov, A. V., Brierley, C. M. & Emanuel, K. A. Tropical cyclones and permanent El Ni ̃no
in the early Pliocene epoch. Nature 463, 1066–1070 (2010).

[9] IPCC. Special Report on the Ocean and Cryosphere in a Changing Climate (SROCC) (September 25, 2019).

[10] Chavas, Daniel, et al. "US hurricanes and economic damage: Extreme value perspective." Natural Hazards Review 14.4 (2013): 237-246.

[11] Marsooli, R., Lin, N., Emanuel, K. and Feng, K., 2019. Climate change exacerbates hurricane flood hazards along US Atlantic and Gulf Coasts in spatially varying patterns. Nature communications, 10(1), pp.1-9.

[12] Marsooli, R. and Lin, N., 2018. Numerical modeling of historical storm tides and waves and their interactions along the US East and Gulf Coasts. Journal of Geophysical Research: Oceans, 123(5), pp.3844-3874.

[13] Lin N, Emmanuel K, 2016. Grey Swan tropical cyclones
