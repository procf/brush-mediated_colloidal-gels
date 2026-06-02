# Brush-mediated angular constraints reshape structure, rigidity, and percolation in colloidal depletion gels
This project is [available on arXiV](https://arxiv.org/abs/2603.13596) and is currently under review.

This is an analysis pipeline for evaluating the effect of non-central, angular constraints in colloidal depletion gels.

Simulations with these angular constraints can be generated using Paniz Haghighi's implementation in HOOMD-blue v4.2.1, available upon request: haghighi.p@northeastern.edu

## What to expect
For systems of ~10,000 colloidal particles, calculate:
- Coordination number distribution and average coordination number
- Void size distribution and average void size
- Network edgelist, the size of connected components, and the physical network diameter ("span" as a proportion of total box size)
- Gaussian Mixture Model (GMM) based mesoscale clustering
- Cauchy-Born estimate of the total elastic modulus from mesoscale (cluster) structure

Note: Standard analyses typically take less than 20min. GMM clustering can take ~9hrs. Cauchy-Born estimate takes <1min.

## Software/package requirements
In this project, the following packages are actively used:
1. GNU Fortran (GCC) 11.4.1
2. `python` v3.10.16 
3. `gsd` v3.2.1
4. `numpy` v1.24.40
5. `pandas` v2.2.3
6. `networkx` v3.4.2
7. `scipy` v1.15.3
8. `node2vec` v0.5.0
9. `umap-learn` v0.5.9
10. `scikit-learn` v1.7.2

## Hardware/OS tested
The program was tested on a single HPC-node running Rocky Linux 9.3 (kernel 5.14).

## Background

Colloidal gel rheology is directly influenced by a gel's underlying particle structure. 
We discovered a new way to modify this structure in experiments. By changing the density of electro-steric surface grafted brush coatings 
on colloidal particles we observe the formation of highly stable, fractal-like gels that exhibit a significantly higher elastic modulus 
(2.7x higher) than their traditional counterparts. This change is consitent with the emergence of non-central angular bending rigidity. 
We developed a new theoretical model that can mimic these effects in simulation. 

The code in this repo quantitatively evaluates the observed change in structure and the mesoscale 
strucutral origins of the change in mechanics.

Here's the abstract: 

*Colloidal gels, like many other soft and disordered solids derive their mechanical properties not only from
the strength of interparticle attraction, but also from the symmetry of the forces that constrain particle motion.
While non-central interactions are known to profoundly alter rigidity and elasticity, they are typically introduced
through particle anisotropy, surface roughness, or patchy interactions, obscuring their independent role. Here we
demonstrate a minimal and geometry-preserving route to emergent non-central forces in colloidal gels by reducing
the density of surface-grafted polymer brushes. At low brush density, partial brush interpenetration introduces
an effective angular bending rigidity at particle contacts, despite fully isotropic particle geometry. This emergent
constraint suppresses local densification, stabilizes low-coordination networks, and produces highly ramified gel
structures with enhanced elasticity. Combining experiments, simulations, and mean-field theory, we show that
these non-central constraints reorganize structure and mechanics across length scales, shifting gelation boundaries
and increasing the elastic modulus by nearly a factor of three. Our results establish surface brush density as a
generic control parameter for programming interaction symmetry in soft particulate matter, with implications for
rigidity, percolation, and mechanical design in disordered systems.*


## Contributors
This is a collaboration with experiments done by the [Soft Matter Engineering Laboratory](https://smel.eng.uci.edu/) at the University of California, Irvine. \
This work was done by Calvin (Ziye) Zhuang, [Rob Campbell](https://scholar.google.com/citations?user=i8S54zYAAAAJ&hl=en), [Paniz Haghighi](https://scholar.google.com/citations?user=LSuqU6YAAAAJ&hl=en), 
[Safa Jamali](https://scholar.google.com/citations?user=D1asaYIAAAAJ&hl=en), and [Ali Mohraz](https://scholar.google.com/citations?user=pW80NaAAAAAJ&hl=en). \
Authors acknowledge support from the National Science Foundation and NASA ROSES FINESST.
