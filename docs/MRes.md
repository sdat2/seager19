# MRes Proposal

No sensitivity analysis is performed in Seager et al. 2019 [1]. Given that the model takes around 30 minutes to run at the moment, it would be quite achievable to do the number of 
model runs neccessary to achieve this. 
Of particular interest might be the model's sensitivity to the 
drag coefficient, as they note that they chose a much 
higher value than normal so as to replicate the amplitude of 
ENSO.
This sensitivity analysis could first be achieved using a Gaussian Process (GP) with a radial basis function (RBF) kernel of a given smoothness,
as the number of data points will initially be quite small (<10^{4}).
From this initial baseline we could expand to more sophisticated
sensitivity analyses, and/or more complicated model settings.

Scientific questions to be addressed would include:
 - Can we replicate the results displayed in Seager et. al 2019?
 - How robust is the model to the parameters chosen?
 - Can the sensitivity of the model to the parameters be understood from the physical processes underlying it?
 - How skillful are different emulation functions at fitting input/output?

Cloud computing platforms like Pangeo could make it relatively easy
to do the same analysis to the CMIP6 ensemble, as they
did to CMIP5 
and see if the bias exists (it seems it does, in e.g. [2]).

Scientific questions to be addressed would include:
 - Do any particular models perform better.
 - What are the likely effects of these biases on risk 
   from e.g. Tropical cyclones (using a metric like the 
   potential intensity).

### Skills learnt, and tie in to AI4ER rubric
This project will be a good opportunity to learn about numerical solutions to differential equations [3, 4]. From the rubric at <https://ai4er-cdt.esc.cam.ac.uk> this project would include the primary application areas of _Weather, Climate and Air Quality_, as well as touching on _Natural Hazards_ as this bias in CMIP will effect the frequency and intensity of tropical cyclones (and therefore their effects such as storm surges, as in my MSci thesis <https://bit.ly/msci-report>).  The methodology of this project will primarily be within _2. Environmental modelling_.  The model used is quite low down in the hierarchy of models, and so it will be much 
easier to demonstrate simple AI tools upon it. 



[1] Seager,  R. et  al.  Strengthening  tropical  Pacific  zonal  sea  surface  temperature  gradient  consistent  withrising  greenhouse  gasesJuly 2019. <https://doi.org/10.1038/s41558-019-0505-x>.

[2] Tian, B. & Dong, X. The Double-ITCZ Bias in CMIP3, CMIP5, and CMIP6 Models Based on Annual Mean Precipitation.Geophysical  Research  Letters 47, e2020GL087232. issn: 0094-8276. doi:10.1029/2020GL087232. <https://onlinelibrary.wiley.com/doi/abs/10.1029/2020GL087232> (Apr. 2020).

[3] Hinch,  E.  J. Think  Before  You  Compute doi:10.1017/9781108855297. <https://www.cambridge.org/core/books/think-before-you-compute/8CC661846951CE3F08F1CC739FBB341F> (Cambridge University Press, Apr. 2020).

[4] Iserles, A. Numerical  Solution  of  Differential  Equations Lent Term 2021 (Centre for Mathematical Sciences, University of Cambridge, 2021)

