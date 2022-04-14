Contains input files to reproduce calculations in [D. Quigley, "Thermodynamics of stacking disorder in ice nuclei", *J. Chem. Phys* 141, 121101 (2014)](https://aip.scitation.org/doi/abs/10.1063/1.4896376). 

In all cases the inputs in `mk` or `mkweights` should be run until a converged `eta_weights.dat` is produced. That file should then be copied to the corresponding `use` or `useweights` folder and the second calculation run until the free energy difference reported in mc.log convergeswith respect to the number of Monte Carlo cycles.
