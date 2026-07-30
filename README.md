# Eigenfunction-expansion-of-an-Oscillating-water-column-slotted-breakwater

[DOI](https://zenodo.org/badge/1317025962.svg)](https://doi.org/10.5281/zenodo.21700609)

Code accompanying the paper

> Reyes, C.C.M. and Huang, Z.
> *An analytical study of a dual-function oscillating water column with nonlinear power take-off and a slotted supporting wall*
> Ocean Engineering (2026)

## Description

This repository contains the numerical implementation of the theory described in the paper: "An analytical study of a dual-function oscillating water column with nonlinear power take-off and a slotted supporting wall" 10.1016/j.oceaneng.2026.127222. 

A non-linear PTO is performed with an iteration procedure performed until convergence of the dissipation parameter \beta_{PTO}. 

Solves the linear dispersion relation to obtain the real and imaginary wave numbers in disp_rel_ee2.m. 
Note that initial guess for the solution must be close to the right hand asymptote (closer to the positive x) and not the left hand. 

---

## Installation

Clone the repository

```bash
git clone https://github.com/clintreyes2/Eigenfunction-expansion-of-an-Oscillating-water-column-slotted-breakwater.git
```

## Citation

If you use this software, please cite

Reyes, C.C.M. 

Eigenfunction expansion of a dual-function oscillating water column-slotted breakwater.

Zenodo

DOI: [https://doi.org/10.5281/zenodo.21700609](https://doi.org/10.5281/zenodo.21700609)

and the accompanying [Ocean Engineering paper](https://www-sciencedirect-com.eres.library.manoa.hawaii.edu/science/article/pii/S0029801826030568?via%3Dihub).

---

## License

MIT License
