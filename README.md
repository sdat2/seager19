# seager19

A repository to contain and analyse the code from:

## Seager et al. 2019, Nature Climate Change

## Strengthening Tropical Pacific Zonal Sea Surface Temperature Gradient Consistent with Rising Greenhouse Gases

### Link:

<https://doi.org/10.1038/s41558-019-0505-x>


### Short summary:


- There is a west-east warm-to-cold contrast in the Pacific.

- GCMs predict weakening contrast with GHG conc.

- In observations it has increased.

- Their cleverer simple linear model agrees with observations.




The citation for this article is:

@article{seager2019strengthening,
  title={Strengthening tropical Pacific zonal sea surface temperature gradient consistent with rising greenhouse gases},
  author={Seager, Richard and Cane, Mark and Henderson, Naomi and Lee, Dong-Eun and Abernathey, Ryan and Zhang, Honghai},
  journal={Nature Climate Change},
  volume={9},
  number={7},
  pages={517--522},
  year={2019},
  url={https://doi.org/10.1038/s41558-019-0505-x},
  publisher={Nature Publishing Group}
}


The paper is discussed at:

<https://deep-convection.org/2020/04/13/episode-5-richard-seager/>



The code was taken from:

<http://kage.ldeo.columbia.edu:81/SOURCES/.LDEO/.ClimateGroup/.PROJECTS/.PublicationsData/.Seager_etal_NCC-2019/>

From running the command

      cloc $(git ls-files)


The initial state off the code was:

      48 text files.
      45 unique files.                              
      14 files ignored.

github.com/AlDanial/cloc v 1.84  T=0.10 s (349.4 files/s, 149217.5 lines/s)

 | Language                |       files       |     blank      |    comment      |       code | 
 | ----------------------- | ----------------- | -------------- | --------------- | ---------- | 
 | Fortran 77              |          15       |      1364      |       1365      |       6170 | 
 | C                       |           5       |       493      |        200      |       2746 | 
 | Jupyter Notebook        |           2       |         0      |        517      |        474 | 
 | Python                  |           2       |       172      |        100      |        397 | 
 | C/C++ Header            |           8       |        88      |         18      |        365 | 
 | make                    |           1       |        15      |          1      |         36 | 
 | Markdown                |           1       |         0      |          0      |          1 | 
 | SUM:                    |          34       |      2132      |       2201      |      10189 | 


## The Beta-Plane Shallow Water Equations

Taken from <https://github.com/milankl/ShallowWaters.jl/edit/master/README.md>


The non-linear shallow water model plus tracer equation is

          ∂u/∂t + (u⃗⋅∇)u - f*v = -g*∂η/∂x - c_D*|u⃗|*u + ∇⋅ν*∇(∇²u) + Fx(x,y)     (1)
          ∂v/∂t + (u⃗⋅∇)v + f*u = -g*∂η/∂y - c_D*|u⃗|*v + ∇⋅ν*∇(∇²v) + Fy(x,y)     (2)
          ∂η/∂t = -∇⋅(u⃗h) + γ*(η_ref - η) + Fηt(t)*Fη(x,y)                       (3)
          ∂ϕ/∂t = -u⃗⋅∇ϕ                                                          (4)

with the prognostic variables velocity u⃗ = (u,v) and sea surface heigth η. The layer thickness is h = η + H(x,y). 
The Coriolis parameter is f = f₀ + βy with beta-plane approximation. 
The graviational acceleration is g. Bottom friction is either quadratic with drag coefficient c_D or linear with inverse time scale r. 
Diffusion is realized with a biharmonic diffusion operator, with either a constant viscosity coefficient ν, 
or a Smagorinsky-like coefficient that scales as ν = c_Smag*|D|, with deformation rate |D| = √((∂u/∂x - ∂v/∂y)² + (∂u/∂y + ∂v/∂x)²). 
Wind forcing Fx is constant in time, but may vary in space.

The linear shallow water model equivalent is

          ∂u/∂t - f*v = -g*∂η/∂x - r*u + ∇⋅ν*∇(∇²u) + Fx(x,y)     (1)
          ∂v/∂t + f*u = -g*∂η/∂y - r*v + ∇⋅ν*∇(∇²v) + Fy(x,y)     (2)
          ∂η/∂t = -H*∇⋅u⃗ + γ*(η_ref - η) + Fηt(t)*Fη(x,y)         (3)
          ∂ϕ/∂t = -u⃗⋅∇ϕ                                           (4)

ShallowWaters.jl discretises the equation on an equi-distant Arakawa C-grid, with 2nd order finite-difference operators. 
Boundary conditions are implemented via a ghost-point copy and each variable has a halo of variable size to account for 
different stencil sizes of various operators.

# Model from the Methods Appendix to Seager 19

## Atmospheric Model

## Ocean Model




