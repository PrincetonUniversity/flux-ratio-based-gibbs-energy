# flux-ratio-based-gibbs-energy
Determine Gibbs free energy of reaction and its confidence interval via metabolic flux analysis and optimization
Integrate absolute cellular metabolite concentrations with flux-derived Gibbs energies and refine confidence intervals

Code used in:
 Conservation of cellular metabolite concentrations and free energies
  Junyoung O. Park, Sara A. Rubin, Yi-Fan Xu, Daniel Amador-Noguez, Jing Fan, Tomer Shlomi, Joshua D. Rabinowitz
   Princeton University

In each organism folder,
- *_mea.xlsx has labeling fractions, measured fluxes, and net/exchange flux inequality constraints that went into the model as input
- .xml shows the metabolites, reactions, carbon mapping, and flux equality constraints
- .h5 and .m function contain the stoichiometric matrix, its kernel, and the cumomer model
- .mat file contains the genome scale model with measured metabolie concentrations and reaction free energies

In brenda folder,
- python script to obtain Km and Ki values from BRENDA

In integrate folder,
- matlab script to solve the quadratic programming (optimal concentration and free energy set) and the linear programming (their refined lower and upper bounds) problems

To run MFA and obtain Gibbs energy of reaction,
1) Open Matlab 2013b or newer
2) Set Path -> Add Folder -> Choose ./src -> Save
3) Open (*_)script.m and execute
4) To run it parallel, enter 'matlabpool local 4' on the command line. '4' can be a different number if a different number of cores is desired
5) Running on clusters requires an additional setup. Please contact system administrator