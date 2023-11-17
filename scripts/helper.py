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
        'M_ac_e': 2,  # Acetate
        'M_acald_e': 2,  # Acetaldehyde
        'M_akg_e': 5,  # Alpha-Ketoglutarate
        'M_co2_e': 1,  # Carbon Dioxide
        'M_etoh_e': 2,  # Ethanol
        'M_for_e': 1,  # Formate
        'M_fru_e': 6,  # Fructose
        'M_fum_e': 4,  # Fumarate
        'M_glc__D_e': 6,  # Glucose
        'M_gln__L_e': 5,  # Glutamine
        'M_glu__L_e': 5,  # Glutamate
        'M_h2o_e': 0,  # Water
        'M_lac__D_e': 3,  # Lactate
        'M_mal__L_e': 4,  # Malate
        'M_nh4_e': 0,  # Ammonium
        'M_o2_e': 0,  # Oxygen
        'M_pi_e': 0,  # Phosphate
        'M_pyr_e': 3,  # Pyruvate
        'M_succ_e': 4  # Succinate
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


def ex_metabolite_formula_lookup():
    return {
        'EX_ac_e': 'CH3COOH',
        'EX_acald_e': 'C2H4O',
        'EX_akg_e': 'C5H6O5',
        'EX_co2_e': 'CO2',
        'EX_etoh_e': 'C2H6O',
        'EX_for_e': 'CH2O2',
        'EX_fru_e': 'C6H12O6',
        'EX_fum_e': 'C4H4O4',
        'EX_glc__D_e': 'C6H12O6',
        'EX_gln__L_e': 'C5H10N2O3',
        'EX_glu__L_e': 'C5H9NO4',
        'EX_h2o_e': 'H2O',
        'EX_lac__D_e': 'C3H6O3',
        'EX_mal__L_e': 'C4H6O5',
        'EX_nh4_e': 'NH4',
        'EX_o2_e': 'O2',
        'EX_pi_e': 'PO₄³⁻',
        'EX_pyr_e': 'C3H4O3',
        'EX_succ_e': 'C4H6O4',
        'EX_h_e': 'H+'
    }


def metabolite_formula_lookup():
    return {
        'M_ac_e': 'CH3COOH',
        'M_acald_e': 'C2H4O',
        'M_akg_e': 'C5H6O5',
        'M_co2_e': 'CO2',
        'M_etoh_e': 'C2H6O',
        'M_for_e': 'CH2O2',
        'M_fru_e': 'C6H12O6',
        'M_fum_e': 'C4H4O4',
        'M_glc__D_e': 'C6H12O6',
        'M_gln__L_e': 'C5H10N2O3',
        'M_glu__L_e': 'C5H9NO4',
        'M_h2o_e': 'H2O',
        'M_lac__D_e': 'C3H6O3',
        'M_mal__L_e': 'C4H6O5',
        'M_nh4_e': 'NH4',
        'M_o2_e': 'O2',
        'M_pi_e': 'PO₄³⁻',
        'M_pyr_e': 'C3H4O3',
        'M_succ_e': 'C4H6O4',
        'M_h_e': 'H+'
    }


def metabolite_to_reaction_name(metabolite_name):
    """
       Convert a metabolite name from a table to a reaction name in a metabolic model.
    """
    if metabolite_name == "objective":
        return "Biomass_Ecoli_core"
    else:
        return "EX_" + metabolite_name.split('_')[1] + "_e"


def carbon_count_for_metabolite(name, val):
    carbon_count_lookup = get_carbon_count_lookup()
    multiplier = carbon_count_lookup.get(name, 0)
    return abs(val) * multiplier


def create_reaction_equations(df: pd.DataFrame, indices: list):
    formulas = []
    lookup = metabolite_formula_lookup()
    for i in indices:
        negativeValues = []
        positiveValues = []
        for key, value in df.loc[i].items():
            if key == 'M_h_e':
                continue

            count_str = ''

            if value < 0:
                value = round(value, 1)
                count_str += str(abs(value)) + ' ' + lookup[key]
                negativeValues.append(count_str)
            elif value > 0:
                value = round(value, 1)
                count_str += str(abs(value)) + ' ' + lookup[key]
                positiveValues.append(count_str)

            if key == 'M_succ_e':
                break  # break after Succinate -> following are calculated rows

        formulas.append('Reaction ' + str(i) + ': ' + ' + '.join(negativeValues) + ' -> ' + ' + '.join(positiveValues))
    return formulas


def create_reaction_equations_atp_biomass(df: pd.DataFrame, df_key: str, indices: list):
    formulas = []
    lookup = ex_metabolite_formula_lookup()
    for i in indices:
        negativeValues = []
        positiveValues = []
        for_dict = df.loc[i, df_key]
        for key, value in for_dict.items():
            if key == 'EX_h_e':
                continue

            count_str = ''

            if value < 0:
                value = round(value, 1)
                count_str += str(abs(value)) + ' ' + lookup[key]
                negativeValues.append(count_str)
            elif value > 0:
                value = round(value, 1)
                count_str += str(abs(value)) + ' ' + lookup[key]
                positiveValues.append(count_str)

        formulas.append('Reaction ' + str(i) + ': ' + ' + '.join(negativeValues) + ' -> ' + ' + '.join(positiveValues))
    return formulas
