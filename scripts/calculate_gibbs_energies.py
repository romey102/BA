import pandas as pd
from scripts.helper import get_atom_sum_for_df


def calculate_gibbs_energies(enthalpies_of_formation_df: pd.DataFrame, results_df: pd.DataFrame) -> pd.DataFrame:
    # Filter out catabolic reactions:
    results_df = results_df[results_df["objective"] == 0]

    # Removing of the column "Protons" & "Biomass":
    results_df = results_df.drop('M_h_e', axis=1)
    results_df = results_df.drop('objective', axis=1)

    # Reset the index:
    enthalpies_of_formation_df = enthalpies_of_formation_df.reset_index(drop=True)
    results_df = results_df.reset_index(drop=True)

    sum_list = []

    for idx, row in results_df.iterrows():
        total_sum = 0
        for col, val in row.items():
            if col in enthalpies_of_formation_df.columns:
                # Directly use the value from the first row in enthalpies_of_formation_df
                if val != 0:
                    corresponding_val = enthalpies_of_formation_df.loc[0, col]
                    total_sum += val * corresponding_val
        sum_list.append(total_sum)
    results_df['Gibbs_Energy'] = sum_list
    return results_df


def standardize_gibbs_energies(gibbs_energy_df: pd.DataFrame) -> pd.DataFrame:
    # Adding the C atom sums as a new column to the DataFrame:
    c_atom_sum_list = get_atom_sum_for_df(gibbs_energy_df)
    gibbs_energy_df['c_atom_sum'] = c_atom_sum_list

    # Create the new column where the "calculated_sum" values are divided by the c_atom_sum_list values:
    gibbs_energy_df['standardized_Gibbs_Energy'] = gibbs_energy_df['Gibbs_Energy'] / c_atom_sum_list
    return gibbs_energy_df
