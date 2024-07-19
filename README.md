### Parametrized quasinormal modes beyond Schwarzschild and Kerr

This repository provides a python wrapper to utilize the parametrized QNM framework. 



#### Schwarzschild black holes

The framework describes linear and quadratic QNM shifts due to modifications of the Regge-Wheeler, Zerilli and scalar perturbations around Schwarzschild black holes beyond general relativity. Currently including complex coefficients for $\alpha_k$ with $k \in [0,\dots, 10]$. 
The following QNMs can be computed:

 - $\ell \in [0,1,2,3,4]$ for scalar and $\ell \in [2,3,4]$ for tensor
 - $n \in [0,1,2]$.



#### Kerr black holes

The framework describes linear QNM shifts and shifts in the separation constant due to modifications of the Teukolsky perturbation equation around Kerr black holes beyond general relativity. Currently including complex coefficients for $\alpha_k$ with $k \in [-6,\dots, 4]$. Note the index shift when specifying the parameter injection that always starts with the zero entry of an array. The following QNMs can be computed:

 - $\ell \in [2,3,4]$
 - $m \in [-\ell, \dots, \ell]$
 - $n \in [0,1,2]$


#### The repository is structured as follows

  - QNM_GR (containts the general relativity values for QNMs)
  - QNM_coefficients (containts the linear/quadratic coefficients)
  - QNM_class.py (defines the QNM class that allows for a quick computation of QNMs)
  - extra_functions.py (some helper functions to load files)
  - tutorial_notebook.ipynb (a quick introduction of how to use the QNM_class)


#### References

The parametrized QNM framework for Schwarzschild has been introduced and extended in:
  - Cardoso et al., Phys.Rev.D 99 (2019) 10, 104077, [arXiv:1901.01265](https://arxiv.org/abs/1901.01265)
  - McManus et al., Phys.Rev.D 100 (2019) 4, 044061, [arXiv:1906.05155](https://arxiv.org/abs/1906.05155)
  - Völkel et al., Phys.Rev.D 105 (2022) 8, 084046, [arXiv:2202.08655](https://arxiv.org/abs/2202.08655)

The parametrized QNM framework for Kerr has been introduced in:
  - Cano et al., in preparation (2024)

The GR values have been taken from the [ringdown website](https://pages.jh.edu/eberti2/ringdown/):
  - Berti et al., Phys.Rev.D 73 (2006) 064030, [arXiv:gr-qc/0512160](https://arxiv.org/abs/gr-qc/0512160)
  - Berti et al., Class.Quant.Grav. 26 (2009) 163001, [arXiv:0905.2975](https://arxiv.org/abs/0905.2975)
