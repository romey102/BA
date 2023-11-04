import os

import pandas as pd
from cobra.io import read_sbml_model

from scripts.equilibrator import calculate_enthalpies_of_formation
from scripts.calculate_gibbs_energies import calculate_gibbs_energies, standardize_gibbs_energies
from scripts.delta_H_delta_S_values import calculate_enthalpy, normalize_gibbs_energy_enthalpy_entropy, \
    visualize_enthalpy, visualize_entropy, visualize_gibbs_energy, calculate_entropy
from scripts.flux_network_reactions import visualize_flux_network_reactions
from cobra_fluxes import get_original_bounds, fba_different_glucose_values, visualize_biomass_vs_glucose, \
    fba_different_oxygen_values, visualize_biomass_vs_oxygen, add_ATP_hydrolysis_reaction, update_exchange_fluxes, \
    restrict_glucose_flow

# setup
os.makedirs('csvs', exist_ok=True)
os.makedirs('intermediate_results', exist_ok=True)
os.makedirs('depictions', exist_ok=True)
core_conversions_df = pd.read_csv('csvs/core_conversions.csv')

# processing

print('calculate_enthalpies_of_formation...')
enthalpies_of_formation_df = calculate_enthalpies_of_formation(core_conversions_df.copy())
enthalpies_of_formation_df.to_csv("csvs/enthalpies_of_formation.csv", index=False)

print('calculate_gibbs_energies...')
gibbs_energy_df = calculate_gibbs_energies(enthalpies_of_formation_df.copy(), core_conversions_df.copy())
gibbs_energy_df.to_csv("intermediate_results/gibbs_energy.csv", index=False)

print('standardize_gibbs_energies...')
standardized_gibbs_energy_df = standardize_gibbs_energies(gibbs_energy_df.copy())
standardized_gibbs_energy_df.to_csv("intermediate_results/standardized_gibbs_energy.csv", index=False)

print('calculate_enthalpy...')
enthalpy_df = calculate_enthalpy(gibbs_energy_df.copy())
enthalpy_df.to_csv("csvs/enthalpy.csv", index=False)

print('calculate_entropy...')
entropy_df = calculate_entropy(enthalpy_df.copy())
entropy_df.to_csv("csvs/gibbs_energy_enthalpy_entropy.csv", index=False)

print('normalize_gibbs_energy_enthalpy_entropy...')
gibbs_energy_enthalpy_entropy_normed_df = normalize_gibbs_energy_enthalpy_entropy(entropy_df.copy())
gibbs_energy_enthalpy_entropy_normed_df.to_csv("csvs/gibbs_energy_enthalpy_entropy_normed.csv", index=False)

print('load model and set original bounds...')
model = read_sbml_model("models/e_coli_core.xml")
original_bounds = get_original_bounds(model)

print('fba_results_for_various_glucose_values...')
fba_different_glucose_values_df = fba_different_glucose_values(model, original_bounds)
fba_different_glucose_values_df.to_csv("csvs/fba_results_for_various_glucose_values.csv", index=False)

print('fba_results_for_various_oxygen_values...')
fba_different_oxygen_values_df = fba_different_oxygen_values(model, original_bounds)
fba_different_oxygen_values_df.to_csv("csvs/fba_results_for_various_oxygen_values.csv", index=False)

print('add ATP hydrolysis as a reaction to model...')
model_with_reaction = add_ATP_hydrolysis_reaction(model)

print('set all exchange fluxes to 0 where metabolite in the ecm tool is 0...')
model_with_reaction_zeroed = update_exchange_fluxes(model_with_reaction)

print('atp_optimal_flow_glucose...')
atp_results_df = restrict_glucose_flow(model_with_reaction_zeroed)
atp_results_df.to_csv("csvs/atp_optimal_flow_glucose.csv", index=False)

# visualizing
print('visualizing...')

visualize_gibbs_energy(gibbs_energy_enthalpy_entropy_normed_df.copy())

visualize_entropy(gibbs_energy_enthalpy_entropy_normed_df.copy())

visualize_enthalpy(gibbs_energy_enthalpy_entropy_normed_df.copy())

visualize_flux_network_reactions(core_conversions_df.copy())

visualize_biomass_vs_glucose(fba_different_glucose_values_df.copy())

visualize_biomass_vs_oxygen(fba_different_oxygen_values_df.copy())
