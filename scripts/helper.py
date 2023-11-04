import pandas as pd


def get_atom_sum_for_df(df: pd.DataFrame) -> list:
    # Lookup of the carbon count:
    carbon_count_lookup = get_carbon_count_lookup()

    # Create an empty list to store the calculated totals for each row:
    c_atom_sum_list = []

    # Traverses each row in the DataFrame:
    for idx, row in df.iterrows():
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
    return c_atom_sum_list

def get_carbon_count_lookup():
    return {
        'M_ac_e': 2,      # Acetate
        'M_acald_e': 2,   # Acetaldehyde
        'M_akg_e': 5,     # Alpha-Ketoglutarate
        'M_co2_e': 1,     # Carbon Dioxide
        'M_etoh_e': 2,    # Ethanol
        'M_for_e': 1,     # Formate
        'M_fru_e': 6,     # Fructose
        'M_fum_e': 4,     # Fumarate
        'M_glc__D_e': 6,  # Glucose
        'M_gln__L_e': 5,  # Glutamine
        'M_glu__L_e': 5,  # Glutamate
        'M_h2o_e': 0,     # Water
        'M_lac__D_e': 3,  # Lactate
        'M_mal__L_e': 4,  # Malate
        'M_nh4_e': 0,     # Ammonium
        'M_o2_e': 0,      # Oxygen
        'M_pi_e': 0,      # Phosphate
        'M_pyr_e': 3,     # Pyruvate
        'M_succ_e': 4     # Succinate
    }


def metabolite_name_lookup():
    return {
        'M_ac_e': 'Acetate',
        'M_acald_e': 'Acetaldehyde',
        'M_akg_e': 'Alpha-Ketoglutarate',
        'M_co2_e': 'Carbon Dioxide',
        'M_etoh_e': 'Ethanol',
        'M_for_e': 'Formate',
        'M_fru_e': 'Fructose',
        'M_fum_e': 'Fumarate',
        'M_glc__D_e': 'Glucose',
        'M_gln__L_e': 'Glutamine',
        'M_glu__L_e': 'Glutamate',
        'M_h2o_e': 'Water',
        'M_lac__D_e': 'Lactate',
        'M_mal__L_e': 'Malate',
        'M_nh4_e': 'Ammonium',
        'M_o2_e': 'Oxygen',
        'M_pi_e': 'Phosphate',
        'M_pyr_e': 'Pyruvate',
        'M_succ_e': 'Succinate',
        'M_h_e': 'H'
    }
