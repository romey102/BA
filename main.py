import os

import pandas as pd

from scripts.equilibrator import calculate_enthalpies_of_formation
from scripts.calculate_gibbs_energies import calculate_gibbs_energies, standardize_gibbs_energies
from scripts.delta_H_delta_S_values import calculate_enthalpy, normalize_gibbs_energy_enthalpy_entropy, \
    visualize_enthalpy, visualize_entropy, visualize_gibbs_energy, calculate_entropy
from scripts.flux_network_reactions import visualize_flux_network_reactions

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

# visualizing
print('visualizing...')

visualize_gibbs_energy(gibbs_energy_enthalpy_entropy_normed_df.copy())

visualize_entropy(gibbs_energy_enthalpy_entropy_normed_df.copy())

visualize_enthalpy(gibbs_energy_enthalpy_entropy_normed_df.copy())

visualize_flux_network_reactions(core_conversions_df.copy())
