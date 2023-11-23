import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scripts.helper import get_carbon_count_lookup, create_reaction_equations


def calculate_enthalpy(core_conversions_g_df: pd.DataFrame) -> pd.DataFrame:
    # Enthalpy lookup:
    delta_H_lookup = {
        'M_ac_e': -486,  # Acetate
        'M_acald_e': -192.2,  # Acetaldehyde (from internet)
        'M_akg_e': -1002.53,  # Alpha-Ketoglutarate (from internet)
        'M_co2_e': -394.1,  # Carbon Dioxide (from book)
        'M_etoh_e': -288.0,  # Ethanol (from book)
        'M_for_e': -410,  # Formate (from book)
        'M_glc__D_e': -1264,  # Glucose (from book)
        'M_glu__L_e': -1003.3,  # Glutamate (from internet)
        'M_h2o_e': -286,  # Water (from book)
        'M_lac__D_e': -687,  # Lactate (from book)
        'M_nh4_e': -133,  # Ammonium (from book)
        'M_o2_e': 0,  # Oxygen (from book)
        'M_pyr_e': -596,  # Pyruvate (from book)
        'M_succ_e': -909  # Succinate (from book)
    }

    # core_conversions_g_df = pd.read_csv("intermediate_results/gibbs_energy.csv", header=0)
    core_conversions_g_df = core_conversions_g_df.reset_index(drop=True)

    # Calculation of Enthalpy:
    delta_H_list = []

    for idx, row in core_conversions_g_df.iterrows():
        total_sum = 0
        for col, val in row.items():
            if col == 'Gibbs_Energy':
                continue
            if col in delta_H_lookup:
                # Directly use the value from the first row in delta_H_df:
                if val != 0:
                    corresponding_val = delta_H_lookup[col]
                    total_sum += val * corresponding_val
        delta_H_list.append(total_sum)

    core_conversions_g_df['Enthalpy'] = delta_H_list
    return core_conversions_g_df


def calculate_entropy(enthalpy_df: pd.DataFrame) -> pd.DataFrame:
    delta_H_G = enthalpy_df
    # Loading the combined data:
    # enthalpy_df = pd.read_csv("csvs/enthalpy_entropy.csv", header=0)

    # Initialize the list for Entropy values:
    delta_S_list = []

    for i in range(len(delta_H_G)):
        delta_H = delta_H_G['Enthalpy'].iloc[i] * 1000  # Conversion from kJ/mol to J/mol
        delta_G = delta_H_G['Gibbs_Energy'].iloc[i] * 1000  # Conversion from kJ/mol to J/mol
        T = 298.15  # Standard temperature in Kelvin
        delta_S = (delta_H - delta_G) / T
        delta_S_list.append(delta_S)

    # Adding the calculated Entropy values to the DataFrame:
    delta_H_G['Entropy'] = delta_S_list

    return delta_H_G


def normalize_gibbs_energy_enthalpy_entropy(gibbs_energy_enthalpy_entropy: pd.DataFrame) -> pd.DataFrame:
    # Normalize all values to C atoms:

    carbon_count_lookup = get_carbon_count_lookup()
    #  Loading the most recent table:
    # delta_H_G_and_S_df = pd.read_csv("csvs/gibbs_energy_enthalpy_entropy.csv", header=0)

    # Create an empty list to store the calculated totals for each row:
    c_atom_sum_list = []

    # Traverses each row in the DataFrame:
    for idx, row in gibbs_energy_enthalpy_entropy.iterrows():
        c_atom_sum = 0  # Initialize the C atom sum for the current row

        # Cycle through each column in the current row:
        for col, val in row.items():
            # Use negative values as positive ones:
            abs_val = abs(val)

            # Using the lookup to find the number of C atoms for the current column:
            if col in carbon_count_lookup:  # Check if the column name exists in the lookup dictionary
                c_atom_sum += carbon_count_lookup[col] * abs_val

        # Add the calculated C atom sum to the list:
        c_atom_sum_list.append(c_atom_sum / 2)

    # Adding the C atom sums as a new column to the DataFrame:
    gibbs_energy_enthalpy_entropy['c_atom_sum'] = c_atom_sum_list

    # Create the new column where the "Gibbs_Energy" values are divided by the c_atom_sum_list values:
    gibbs_energy_enthalpy_entropy['Gibbs_Energy'] = gibbs_energy_enthalpy_entropy[
                                                        'Gibbs_Energy'] / c_atom_sum_list  # kJ/mol*C
    gibbs_energy_enthalpy_entropy['Enthalpy'] = gibbs_energy_enthalpy_entropy['Enthalpy'] / c_atom_sum_list  # kJ/mol*C
    gibbs_energy_enthalpy_entropy['Entropy'] = gibbs_energy_enthalpy_entropy['Entropy'] / c_atom_sum_list  # J/K*C

    # Save the updated DataFrame as a CSV file:
    return gibbs_energy_enthalpy_entropy


def visualize_gibbs_energy(gibbs_energy_enthalpy_entropy_normed: pd.DataFrame) -> None:
    formulas = create_reaction_equations(
        gibbs_energy_enthalpy_entropy_normed.sort_values('Gibbs_Energy', ascending=True, ignore_index=True),
        [0, 1, 2, 3, 4, 5])

    legend_elements = [Line2D([0], [0], color='w', marker='o', markerfacecolor='w', markersize=10, label=sf) for sf in
                       formulas]

    gibbs_energy_enthalpy_entropy_normed = gibbs_energy_enthalpy_entropy_normed.sort_values(by='Gibbs_Energy',
                                                                                            ascending=True)

    plt.figure(figsize=(14, 7))
    barWidth = 0.3

    # Use of absolute values for Gibbs Energy:
    plt.bar(range(len(gibbs_energy_enthalpy_entropy_normed)),
            gibbs_energy_enthalpy_entropy_normed['Gibbs_Energy'].abs(), color='r', label='Gibbs Energy')

    # Adding labels for the bars:
    plt.xlabel('Reactions', fontweight='bold')
    plt.ylabel('ΔG [- kJ/C-mol]')
    plt.title('Standardized Gibbs energy changes for the catabolic reactions of the E. coli core model')
    plt.xticks(range(len(gibbs_energy_enthalpy_entropy_normed)), range(len(gibbs_energy_enthalpy_entropy_normed)))

    plt.legend(handles=legend_elements, loc="upper right")

    # Save diagram:
    print(f"Saving gibbs_energy.png")
    plt.savefig('depictions/gibbs_energy.png')
    plt.show()


def visualize_entropy(gibbs_energy_enthalpy_entropy_normed: pd.DataFrame) -> None:
    formulas = create_reaction_equations(
        gibbs_energy_enthalpy_entropy_normed.sort_values('Gibbs_Energy', ascending=True, ignore_index=True),
        [0, 1, 2, 3, 4, 5])

    legend_elements = [Line2D([0], [0], color='w', marker='o', markerfacecolor='w', markersize=10, label=sf) for sf in
                       formulas]

    gibbs_energy_enthalpy_entropy_normed = gibbs_energy_enthalpy_entropy_normed.sort_values(by='Gibbs_Energy',
                                                                                            ascending=True)

    average_of_entropy = gibbs_energy_enthalpy_entropy_normed['Entropy'].mean()
    # Create a color list based on the values of Entropy:
    colors = ['r' if value < 0 else 'g' for value in gibbs_energy_enthalpy_entropy_normed['Entropy']]

    plt.figure(figsize=(14, 7))
    barWidth = 0.3

    # Using the color list when drawing the bars:
    plt.bar(range(len(gibbs_energy_enthalpy_entropy_normed)), gibbs_energy_enthalpy_entropy_normed['Entropy'],
            color=colors, label='Entropy')

    plt.axhline(y=average_of_entropy, color='b', linestyle='-', label=f'Sum of Entropy: {average_of_entropy}')

    # Adding labels for the bars:
    plt.xlabel('Reactions', fontweight='bold')
    plt.ylabel('ΔS [J/T*C]')
    plt.title('Standardized entropy changes for the catabolic reactions of the E. coli core model')
    plt.xticks(range(len(gibbs_energy_enthalpy_entropy_normed)), range(len(gibbs_energy_enthalpy_entropy_normed)))

    plt.legend(handles=legend_elements + [
        Line2D([0], [0], color='b', linestyle='-', label=f'Sum of Entropy: {average_of_entropy}')], loc="upper right")

    plt.legend(handles=legend_elements, loc="upper right")

    # Save diagram:
    print(f"Saving entropy.png")
    plt.savefig('depictions/entropy.png')
    plt.show()


def visualize_enthalpy(gibbs_energy_enthalpy_entropy_normed: pd.DataFrame) -> None:
    formulas = create_reaction_equations(
        gibbs_energy_enthalpy_entropy_normed.sort_values('Gibbs_Energy', ascending=True, ignore_index=True),
        [0, 1, 2, 3, 4, 5])

    legend_elements = [Line2D([0], [0], color='w', marker='o', markerfacecolor='w', markersize=10, label=sf) for sf in
                       formulas]

    gibbs_energy_enthalpy_entropy_normed = gibbs_energy_enthalpy_entropy_normed.sort_values(by='Gibbs_Energy',
                                                                                            ascending=True)

    plt.figure(figsize=(14, 7))
    barWidth = 0.3

    # Use of absolute values for Enthalpy:
    plt.bar(range(len(gibbs_energy_enthalpy_entropy_normed)), gibbs_energy_enthalpy_entropy_normed['Enthalpy'].abs(),
            color='r', label='Enthalpy')

    # Adding labels for the bars:
    plt.xlabel('Reactions', fontweight='bold')
    plt.ylabel('ΔH [- kJ/C-mol]')
    plt.title('Standardized enthalpy changes for the catabolic reactions of the E. coli core model')
    plt.xticks(range(len(gibbs_energy_enthalpy_entropy_normed)), range(len(gibbs_energy_enthalpy_entropy_normed)))

    plt.legend(handles=legend_elements, loc="upper right")

    # Save diagram:
    print(f"Saving enthalpy.png")
    plt.savefig('depictions/enthalpy.png')
    plt.show()
