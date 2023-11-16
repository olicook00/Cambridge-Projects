# Modelling-a-Hydrogenator
In this project I modelled a hydrogenator, which takes a mix of eight unsaturated hydrocarbons and saturates them. The reactor is a fixed bed reactor and was modelled as an ideal plug flow reactor. The incoming reactant concentrations are fixed from upstream and the reaction kinetics were taken from literature, from which a system of PDE's was formed. Since hydrogen is involved in all the reactions, the PDE's must be solved using a numerical method. Time is discritesed into small steps and concentration as a function of residence time is found using the finite difference method.

Concentration of the reactants as a function of residence time is then plotted. In particular, the reactant 1,3-butadiene had to be reduced to acceptable levels, and this concentration is also shown on the plot.

Assumptions made in the model: 

-  Ideal plug flow reactor. Unlikely, but much simpler than producing a residence time distribution.
- The total molar concentration is constant. Therefore, molar concentration of each reactant is only a function of its mol fraction, which simplifies the system significantly. Actually accurate to ~2% according to UniSim thermodynamics.
-  The thermodynamic behaviour was very simplified: hydrogen is assumed to be completely soluble at this pressure and concentration
