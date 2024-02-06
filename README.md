# “Systematic computational characterization of all possible catabolic routes of the E.coli core model"
Bachelor thesis by Ronja Meyer

---

In this Python project for my bachelor thesis, I am taking a closer look at the catabolic metabolic pathways of E. coli. 
For this, I have written various calculations that can be used to compute and plot the enthalpies of formation, Gibbs energies, enthalpies, entropies and their standardized values, 
as well as using CobraPy to calculate, standardize and plot the FBAs for biomass at different glucose and oxygen bounds, calculate the optimal ATP flux for glucose, and compute, standardize and plot the maximum ATP and biomass yield for the reactions. 
This calculations require results from the ecmtool for the model (in this case the E coli core model)! 
You can find my results for the individual calculations in the wiki.
There is also some documentation in the code to make it easier to understand.
For a detailed explanation of the procedure, please refer to my bachelor thesis of the same name in the methods section.

---
# Results:
The different flow diagrams for each individual catabolic reaction:
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_29.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_53.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_55.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_57.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_58.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_59.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_60.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_62.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_131.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_138.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_158.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_160.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_161.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_241.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_242.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_243.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_244.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_245.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_246.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_247.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_248.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_249.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_250.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_251.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_252.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_253.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_254.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_255.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_256.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_257.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_258.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_259.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_260.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_261.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_262.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_263.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_264.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_265.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_266.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_267.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_268.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_269.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_270.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_271.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_272.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_273.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_274.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_275.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_276.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_277.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_461.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_464.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_466.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_469.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_471.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_472.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_473.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_474.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_475.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_485.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_486.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_488.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_489.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_490.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_491.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_492.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_493.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_494.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_495.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_682.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_684.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_685.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_686.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_687.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_688.png "")
![](C:\Users\Ronja\PycharmProjects\BA\depictions\flux_network_reactions/flux_network_reaction_689.png "")

---

The standardized Gibbs energy values for each catabolic reaction:
![](C:\Users\Ronja\PycharmProjects\BA\depictions/gibbs_energy.png "")

The standardized enthalpy values for each catabolic reaction:
![](C:\Users\Ronja\PycharmProjects\BA\depictions/enthalpy.png "")

The standardized entropy values for each catabolic reaction (blue line: average):
![](C:\Users\Ronja\PycharmProjects\BA\depictions/entropy.png "")

---
The standardized biomass production for different standardized glucose lower bounds:
![](C:\Users\Ronja\PycharmProjects\BA\depictions/standardized_biomass_vs_glucose.png "")

The standardized biomass production for different oxygen lower bounds:
![](C:\Users\Ronja\PycharmProjects\BA\depictions/standardized_biomass_vs_oxygen.png "")

---
The standardized maximum ATP values for each catabolic reaction:
![](C:\Users\Ronja\PycharmProjects\BA\depictions/standardized_max_ATP.png "")

The standardized maximum biomass values for each catabolic reaction:
![](C:\Users\Ronja\PycharmProjects\BA\depictions/standardized_max_biomass.png "")
