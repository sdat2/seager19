# seager19

A repository to contain and analyse the code from:

## Seager et al. 2019, Nature Climate Change

## Strengthening Tropical Pacific Zonal Sea Surface Temperature Gradient Consistent with Rising Greenhouse Gases

### Link:

<https://doi.org/10.1038/s41558-019-0505-x>


### Short summary of paper:


- There is a west-east warm-to-cold contrast in the Pacific.

- GCMs predict weakening contrast with GHG conc.

- In observations it has increased.

- Their cleverer simple linear model agrees with observations.


The citation for this paper is:

```
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
```

The paper is discussed in a podcast available at:

<https://deep-convection.org/2020/04/13/episode-5-richard-seager/>


The code and data was taken from:

<http://kage.ldeo.columbia.edu:81/SOURCES/.LDEO/.ClimateGroup/.PROJECTS/.PublicationsData/.Seager_etal_NCC-2019/>

The Python code for the atmosphere model is in a Juypter Notebook. The ocean model code is built on legacy Fortran 90 and C code.

The data is currently not stored in the github repository, as it takes up roughly 3.5 GB.

### Code Makeup

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


### Code structure 

The code is structured into folders:

```
   |-animations
   |-atmos-model
   |---README.md --> lists file structure of this model.
   |---DATA
   |---tmp
   |-ocean-model
   |---README.md --> lists file structure of this model.
   |---DATA
   |---RUN
   |-----run-model --> run this with bash to run model?
   |-----DATA
   |-----output
   |---SRC
   |-----DATA
   |-----output
   |-requirements
```

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


# Model from the Methods Appendix to Seager et al. 2019

## Atmos and Ocean reanalysis

- ECMWF ERA-40 1958-78 & ERA-Interim 1979-2017
  - Wind at 2m from surface, 
  - precipitation, 
- ORAS4 1958-2017
  - SST, 
  - themocline depth (20C isotherm).

## Atmospheric Model

   (u, v, w) = (u'cos(pi z/zT), v'cos(pi z/zT), w'sin(pi z/zT))

z is the height, and zT is the top of the troposphere.

  (theta, Q) = (theta', Q')(theta0/theta00) sin(pi z/zT)

  p = p'(rho0 / rho00) cos(pi z/zT)

  (p/rho0)_z = g theta / theta00

  p'/rho00 = (g zT/ pi theta00)theta'


theta0 and rho0 are basic-state potential temperature and density profiles.

  epsilon_u u - fv + phi_x = 0
  epsilon_v v + fu + phi_y = 0
  epsilon_phi phi + u_x + v_y = -Q_1

the geopotential:

  phi = - (g zT/ pi theta_00) theta

## Ocean Model


## Compilers in Jasmin for fortran:

- `f95`
  ```
  Using built-in specs.
  COLLECT_GCC=f95
  COLLECT_LTO_WRAPPER=/usr/libexec/gcc/x86_64-redhat-linux/4.8.5/lto-wrapper
  Target: x86_64-redhat-linux
  Configured with: ../configure --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-linker-build-id --with-linker-hash-style=gnu --enable-languages=c,c++,objc,obj-c++,java,fortran,ada,go,lto --enable-plugin --enable-initfini-array --disable-libgcj --with-isl=/builddir/build/BUILD/gcc-4.8.5-20150702/obj-x86_64-redhat-linux/isl-install --with-cloog=/builddir/build/BUILD/gcc-4.8.5-20150702/obj-x86_64-redhat-linux/cloog-install --enable-gnu-indirect-function --with-tune=generic --with-arch_32=x86-64 --build=x86_64-redhat-linux
  Thread model: posix
  gcc version 4.8.5 20150623 (Red Hat 4.8.5-39) (GCC) 
  ```

 - `gfortran44`
  ```
  Using built-in specs.
  Target: x86_64-redhat-linux
  Configured with: ../configure --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-languages=c,c++,fortran --disable-libgcj --without-cloog --without-ppl --with-tune=generic --with-arch_32=i686 --build=x86_64-redhat-linux
  Thread model: posix
  gcc version 4.4.7 20120313 (Red Hat 4.4.7-8) (GCC) 
  ```

- `gfortran`
  ```
  Using built-in specs.
  COLLECT_GCC=gfortran
  COLLECT_LTO_WRAPPER=/usr/libexec/gcc/x86_64-redhat-linux/4.8.5/lto-wrapper
  Target: x86_64-redhat-linux
  Configured with: ../configure --prefix=/usr --mandir=/usr/share/man --infodir=/usr/share/info --with-bugurl=http://bugzilla.redhat.com/bugzilla --enable-bootstrap --enable-shared --enable-threads=posix --enable-checking=release --with-system-zlib --enable-__cxa_atexit --disable-libunwind-exceptions --enable-gnu-unique-object --enable-linker-build-id --with-linker-hash-style=gnu --enable-languages=c,c++,objc,obj-c++,java,fortran,ada,go,lto --enable-plugin --enable-initfini-array --disable-libgcj --with-isl=/builddir/build/BUILD/gcc-4.8.5-20150702/obj-x86_64-redhat-linux/isl-install --with-cloog=/builddir/build/BUILD/gcc-4.8.5-20150702/obj-x86_64-redhat-linux/cloog-install --enable-gnu-indirect-function --with-tune=generic --with-arch_32=x86-64 --build=x86_64-redhat-linux
  Thread model: posix
  gcc version 4.8.5 20150623 (Red Hat 4.8.5-39) (GCC) 
  ```

## Current compilation:

```
gfortran  -ffixed-line-length-132 -Ofast -c om_core.F
gfortran  -ffixed-line-length-132 -Ofast -c om_ekm.F
gfortran  -ffixed-line-length-132 -Ofast -c om_sst.F
gfortran  -ffixed-line-length-132 -Ofast -c om_leap.F
gfortran  -ffixed-line-length-132 -Ofast -c om_equi.F
gfortran  -ffixed-line-length-132 -Ofast -c om_forc.F
gfortran  -ffixed-line-length-132 -Ofast -c om_qflux.F
gfortran  -ffixed-line-length-132 -Ofast -c om_tios.F
gfortran  -ffixed-line-length-132 -Ofast -c om_mem.F
gfortran  -ffixed-line-length-132 -Ofast -c om_wrap.F
gfortran  -ffixed-line-length-132 -Ofast -c fodb.F
gfortran  -ffixed-line-length-132 -Ofast -c om_main.F
gcc -Wno-implicit-function-declaration -Ofast -c om_c.c
gcc -Wno-implicit-function-declaration -Ofast -c codb.c
gcc -Wno-implicit-function-declaration -Ofast -c daio.c
gcc -Wno-implicit-function-declaration -Ofast -c sio.c
gfortran  -ffixed-line-length-132 -Ofast  -o tcom om_main.o wrap-mod.o data-mod.o sst-mod.o om_core.o om_ekm.o om_sst.o om_leap.o om_equi.o om_forc.o om_qflux.o om_tios.o om_mem.o om_wrap.o fodb.o  om_c.o codb.o daio.o sio.o -lnetcdf -lnetcdff
```

### Commands to produce file trees:

```
find . -not -path '*/\.*' | python -c "import sys as s;s.a=[];[setattr(s,'a',list(filter(lambda p: c.startswith(p+'/'),s.a)))or (s.stdout.write('  '*len(s.a)+c[len(s.a[-1])+1 if s.a else 0:])or True) and s.a.append(c[:-1]) for c in s.stdin]"
```

```
ls -R | grep ":$" | sed -e 's/:$//' -e 's/[^-][^\/]*\//--/g' -e 's/^/   /' -e 's/-/|/'
```
## Youtube links:

https://doi.org/10.1038/s41558-019-0505-x

```
ocean-model/RUN/output/

SST_SST full in om_run2f: <https://youtu.be/JA97IWPmwxs>
DYN_PRES in om_run2f: <https://youtu.be/5oRMWWAK1sM>
TDEEP_HMODEL in om_run2f: <https://youtu.be/n25l6uYWEzY>
TDEEP_HTHERM in om_run2f: <https://youtu.be/ikOo6VTXfkg>
TDEEP_TDEEP in om_run2f: <https://youtu.be/BSRyTuESzLA>
```