# GExt.jl: Grassmann Extrapolation of Density Matrices for Born–Oppenheimer Molecular Dynamics
Authors: Étienne Polack, Geneviève Dusson, Benjamin Stamm and Filippo Lipparini

Julia template for the Grassmann Extrpolation method for density matrices ([JCTC](https://doi.org/10.1021/acs.jctc.1c00751), [arXiv](https://arxiv.org/abs/2107.13218))

## Abstract

_Born–Oppenheimer molecular dynamics (BOMD) is a powerful but expensive technique. The main bottleneck in a density functional theory BOMD calculation is the solution to the Kohn–Sham (KS) equations that requires an iterative procedure that starts from a guess for the density matrix. Converged densities from previous points in the trajectory can be used to extrapolate a new guess; however, the nonlinear constraint that an idempotent density needs to satisfy makes the direct use of standard linear extrapolation techniques not possible. In this contribution, we introduce a locally bijective map between the manifold where the density is defined and its tangent space so that linear extrapolation can be performed in a vector space while, at the same time, retaining the correct physical properties of the extrapolated density using molecular descriptors. We apply the method to real-life, multiscale, polarizable QM/MM BOMD simulations, showing that sizeable performance gains can be achieved, especially when a tighter convergence to the KS equations is required._

## Citation

```bibtex
@article{doi:10.1021/acs.jctc.1c00751,
  author       = {Polack, Étienne and Dusson, Geneviève and Stamm, Benjamin and Lipparini, Filippo},
  url          = {https://doi.org/10.1021/acs.jctc.1c00751},
  date         = {2021},
  doi          = {10.1021/acs.jctc.1c00751},
  eprint       = {https://doi.org/10.1021/acs.jctc.1c00751},
  journaltitle = {Journal of Chemical Theory and Computation},
  note         = {PMID: 34623810},
  number       = {11},
  pages        = {6965--6973},
  title        = {Grassmann Extrapolation of Density Matrices for Born–Oppenheimer Molecular Dynamics},
  volume       = {17},
}
```
