# Computational Astrophysics and Statistics Projects

This repository contains two course projects completed as part of the MSc course **Computational Astrophysics and Statistics** at the University of Bologna.

The projects provided hands-on training in numerical astrophysics, including finite difference methods, spherical grids, diffusion equations, hydrodynamical simulations, numerical stability, code validation, and the interpretation of numerical results in an astrophysical context.

The repository contains the final reports, numerical codes, and plotting scripts used in the two projects.

## Projects

### Project 1: Hydrostatic Equilibrium and Turbulent Metal Diffusion in Galaxy Clusters

**Report:** [Project1_Wang.pdf](./Project1_Wang.pdf)

This project investigates the structure and chemical evolution of the intracluster medium in a galaxy cluster.

The gas density distribution was obtained by solving the hydrostatic equilibrium equation in a gravitational potential including a Navarro Frenk White dark matter halo and a Hernquist model for the central brightest cluster galaxy. Both isothermal and radially varying gas temperature profiles were considered.

The second part of the project models the evolution of iron abundance in the intracluster medium. Iron transport was described through turbulent diffusion, while enrichment from Type Ia supernovae and stellar winds was included through source terms.

Several physical scenarios were explored, including diffusion only, source only, and combined diffusion plus source models. The resulting metallicity profiles were compared with an observationally motivated Perseus cluster abundance profile.

#### Numerical methods

The main numerical components of this project include:

* One dimensional spherical radial grids
* Staggered placement of cell centered and interface quantities
* Numerical integration of enclosed mass profiles in spherical shells
* Hydrostatic equilibrium calculations in combined gravitational potentials
* Explicit finite difference integration of the metal diffusion equation
* Forward Time Centered Space method
* Stability controlled timesteps for the diffusion calculation
* Outflow boundary conditions
* Stellar wind and Type Ia supernova source terms
* Conservation checks and comparison with analytic solutions
* Parameter studies of the turbulent diffusion coefficient and supernova enrichment rate

The project also included comparisons between numerical and analytic mass and gas density profiles as basic validation tests.

#### Main physical application

The simulations were used to investigate how turbulent transport and stellar injection jointly shape the radial metallicity distribution of the intracluster medium.

The results illustrate the competition between central metal injection and outward turbulent mixing, and show how variations in the diffusion coefficient and supernova enrichment rate modify the resulting abundance profile.


### Project 2: One Dimensional ZEUS Hydrodynamics and Cluster Cooling Flows

**Report:** [Project2_Wang.pdf](./Project2_Wang.pdf)

This project develops a one dimensional hydrodynamical solver based on the numerical framework of the ZEUS code described by **[Stone & Norman (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJS...80..753S)**.

The solver was implemented for both Cartesian and spherical geometries and was first tested using standard hydrodynamical benchmark problems. It was then applied to the evolution of radiatively cooling gas in a galaxy cluster.

The cluster model includes gravitational contributions from a Navarro Frenk White dark matter halo, a Hernquist stellar component, and a central supermassive black hole. Radiative cooling, stellar mass return, Type Ia supernova enrichment, and stellar energy injection were included as source terms.

#### Numerical methods

The hydrodynamics implementation includes:

* One dimensional Eulerian finite difference hydrodynamics
* Staggered computational meshes
* Cartesian and spherical coordinate systems
* Operator splitting between source and transport steps
* Pressure gradient and gravitational source terms
* Von Neumann Richtmyer artificial viscosity for shock treatment
* First order upwind advection
* Numerical transport of mass, internal energy, momentum, and iron abundance
* Courant Friedrichs Lewy timestep control
* Reflective and outflow boundary conditions
* Uniform and nonuniform radial grids
* Radiative cooling source terms
* Stellar mass, energy, and metal injection

The implementation follows the basic ZEUS methodology introduced by [Stone & Norman (1992)](https://ui.adsabs.harvard.edu/abs/1992ApJS...80..753S), adapted here to a one dimensional educational hydrodynamics code.

#### Code validation

Before applying the solver to cluster cooling flows, the numerical implementation was tested using:

* Cartesian Sod shock tube tests
* Spherical shock tube tests
* Strong shock tests
* Resolution comparisons using 100 and 1000 computational zones

The tests reproduce the characteristic rarefaction wave, contact discontinuity, and shock structure expected from standard shock tube solutions.

Increasing the numerical resolution produces sharper discontinuities and improved agreement with reference solutions, providing a basic validation of the hydrodynamical implementation.

#### Cooling flow application

The validated hydrodynamics solver was applied to a one dimensional galaxy cluster cooling flow model.

The simulations follow the evolution of:

* Gas density
* Radial velocity
* Temperature
* Pressure
* Hot and cold gas mass
* Iron abundance
* X-ray surface brightness
* X-ray luminosity
* Mass cooling rate

The model produces inward gas motions driven by radiative cooling and evolves toward a quasi steady cooling flow configuration.

Simulated cooling rates were also compared with analytic estimates derived from the X-ray luminosity.


## Repository Structure

A possible organization of the repository is:

```text
computational-astrophysics-projects/
│
├── README.md
│
├── project1_icm_diffusion/
│   ├── README.md
│   ├── Project1_Wang.pdf
│   ├── src/
│   ├── plotting/
│   └── figures/
│
└── project2_zeus_hydrodynamics/
    ├── README.md
    ├── Project2_Wang.pdf
    ├── src/
    ├── plotting/
    └── figures/
```
The `src/` directories contain the numerical implementations, while `plotting/` contains scripts used to analyze the simulation output and generate the figures shown in the reports.

### Contribution Statement
These projects were completed jointly by Jintong Wang and Jiya Yao.

For both projects, the two authors independently worked through the full set of numerical tasks and subsequently cross checked the implementations and numerical results.

For Report 1, I took primary responsibility for the latter part of Project 1, particularly the turbulent diffusion, source modeling, and numerical analysis. For Report 2, I completed and documented the majority of the project 2, including most of the numerical implementation, hydrodynamical testing, cooling flow analysis, and report preparation. The final results and reports were cross checked by both authors.

The codes included in this repository document my implementation and analysis unless otherwise stated.

## Limitations

These projects were designed as training projects in computational astrophysics rather than as original numerical method research.

The main limitations are:
- Both simulations are one dimensional.
- The ZEUS type solver uses a first order upwind transport method and artificial viscosity rather than a modern high order Godunov scheme.
- No adaptive mesh refinement is implemented.
- No magnetohydrodynamics is included.
- No multidimensional fluid instabilities or turbulence are resolved directly.
- The turbulent metal transport in Project 1 is represented using an effective diffusion coefficient rather than explicit turbulent hydrodynamics.
- The cooling flow model does not include AGN feedback, thermal conduction, or other heating processes required for a fully realistic cluster core.
- The numerical tests provide basic validation and resolution comparisons, but they are not intended as a formal convergence study of the numerical method.

The purpose of the projects was to develop a practical understanding of numerical discretization, hydrodynamical algorithms, stability conditions, validation tests, and astrophysical applications.

## Skills Demonstrated
These projects provided practical experience with:
### Programming
- FORTRAN
- Python for analysis and visualization
  
### Numerical methods
- Finite difference methods
- Explicit time integration
- Staggered grids
- Operator splitting
- Upwind advection
- Artificial viscosity
- CFL stability conditions
- Numerical boundary conditions
- Spherical coordinate discretization

### Computational astrophysics
- Hydrostatic equilibrium
- Galaxy cluster gravitational potentials
- Turbulent diffusion
- Compressible hydrodynamics
- Shock propagation
- Radiative cooling
- Chemical enrichment
- X-ray observables
- Code verification and validation
  
## Notes
These projects were completed as coursework
