# MRes Proposal

![SST output over time period](gifs/SST_SST2_in_om_run2f.gif)

## Part 1

Seager et al. 2019 [1] shows that CMIP models do not follow the observed trend in ENSO, but that this trend can
be reproduced with a simple coupled physical model.
There was no sensitivity analysis of the model in this paper. Given that the model 
takes around 30 minutes to run at the moment, 
it would be quite achievable to do the number of 
model runs neccessary to achieve this. 
Of particular interest might be the model's sensitivity to the 
drag coefficient, as they note that they chose a much 
higher value than normal so as to replicate the amplitude of ENSO.
This sensitivity analysis could first be achieved using a Gaussian Process (GP)
with a radial basis function (RBF) kernel of a given smoothness,
as the number of data points will initially be quite small (<10^{4}).
From this initial baseline we could expand to more sophisticated
sensitivity analyses, and/or more complicated model settings.

Scientific questions to be addressed would include:

 - Can we replicate the results displayed in Seager et al. 2019?
 - How robust is the model to the parameters chosen?
 - Can the sensitivity of the model to the parameters be understood from the physical processes underlying it?
 - How skillful are different emulation functions at fitting input/output of the model?

## Part 2

Cloud computing platforms like Pangeo could make it relatively easy
to do the same analysis to the CMIP6 ensemble, as they
did to CMIP5,
to see if the same bias exists (it seems it does, in e.g. [2]).

Scientific questions to be addressed would include:

 - Do any particular model families perform better, and is there a physical reason for this.
 - What are the likely effects of these biases on risk 
   from e.g. Tropical cyclones (using a metric like the 
   potential intensity).

## Skills learnt

This project will be a good opportunity to learn about numerical solutions to differential equations [3, 4]. From the rubric at <https://ai4er-cdt.esc.cam.ac.uk> this project would include the primary application areas of _Weather, Climate and Air Quality_, as well as touching on _Natural Hazards_ as this bias in CMIP will effect the frequency and intensity of tropical cyclones (and therefore their effects such as storm surges, as in my MSci thesis <https://bit.ly/msci-report>).  The methodology of this project will primarily be within _2. Environmental modelling_.  The model used is quite low down in the hierarchy of models, and so it will be much 
easier to demonstrate simple ML tools upon it. 


### Citations

[1] Seager,  R. et  al.  Strengthening  tropical  Pacific  zonal  sea  surface  temperature  gradient  consistent  withrising  greenhouse  gases July 2019. <https://doi.org/10.1038/s41558-019-0505-x>.

[2] Tian, B. & Dong, X. The Double-ITCZ Bias in CMIP3, CMIP5, and CMIP6 Models Based on Annual Mean Precipitation.Geophysical  Research  Letters 47, e2020GL087232. issn: 0094-8276. doi:10.1029/2020GL087232. <https://onlinelibrary.wiley.com/doi/abs/10.1029/2020GL087232> (Apr. 2020).

[3] Hinch,  E.  J. Think  Before  You  Compute doi:10.1017/9781108855297. <https://www.cambridge.org/core/books/think-before-you-compute/8CC661846951CE3F08F1CC739FBB341F> (Cambridge University Press, Apr. 2020).

[4] Iserles, A. Numerical  Solution  of  Differential  Equations Lent Term 2021 (Centre for Mathematical Sciences, University of Cambridge, 2021)

![Hurricane Katrina at Landfall: Credit NASA](https://cdn.britannica.com/74/121674-050-C458B2B5/satellite-image-National-Oceanic-and-Atmospheric-Administration-August-28-2005.jpg)

## Potential PhD project.

The current deliberately broad title of this project would be:

_Investigating the bias in the predicted risk from Tropical Cyclones_

This could focus on:


(1)  The bias in GCMs, and how this would be expected to effect broad proxies 
   for tropical cyclone activity, such as potential intensity.

   _and / or_
  
 (2) A more focussed look at how the risk from tropical cyclones
   can be parameterised from its characteristics and the characteristics of the coast, 
   involving more extensive extreme value theory.

To achieve this I am in contact with Talea Mayo (storm surge expert), 
and Dan Chavas (tropical cyclone expert).

### Focus 1

The potential intensity index is a good indicator of the maximum size that a storm
could reach in a given climate, which is dependent on the sea surface temperature and
the temperature of the troposphere.
There are also various cyclone genesis indices.
The biases in the global climate models are non stationary.
Correcting the bias in global climate models would require 
some understanding of what the true sentivities of the climate
to the forcing. Perhaps a combination of physical models 
and ML emulators can improve this understanding.

These biases probably arise from errors in the coupled parameterisations,
and errors in the parameterisation of convection. There are probably many
interesting numerical experiments that could be done to investigate this.


### Focus 2

The second focus would be more in line 
with Talea Mayo's work, as she has a number 
of prerun high resolution model outputs 
from tropical cyclones hitting the coast.
My MSci Project <https://bit.ly/msci-report> 
showed that you can predict the responsiveness of
a coastline to a windstess fairly well with a very simple physical model.
To properly cacluate the return period for a given coastline,
quite a lot of data is required. 
In fact, if we consider a single point at the coast, it 
may take roughly 1 million years to properly converge to an extreme value
distribution. Theforefore, ways to pool data from different points along the 
coast would probably be neccessary to be able to achieve anything reliable 
within reasonable computational expense.
