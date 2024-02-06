from multiprocessing import freeze_support

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from cobra import Model
from matplotlib.lines import Line2D
from scripts.helper import carbon_count_for_metabolite, metabolite_to_reaction_name, \
    create_reaction_equations_atp_biomass, get_carbon_count_lookup_biomass, get_carbon_count_lookup


def get_original_bounds(model):
    # List of reactions for which the bounds are stored:
    reaction_list = ['ACALD', 'ACALDt', 'ACKr', 'ACONTa', 'ACONTb', 'ACt2r', 'ADK1', 'AKGDH', 'AKGt2r', 'ALCD2x',
                     'ATPM',
                     'ATPS4r', 'BIOMASS_Ecoli_core_w_GAM', 'CO2t', 'CS', 'CYTBD', 'D_LACt2', 'ENO', 'ETOHt2r',
                     'EX_ac_e',
                     'EX_acald_e', 'EX_akg_e', 'EX_co2_e', 'EX_etoh_e', 'EX_for_e', 'EX_fru_e', 'EX_fum_e',
                     'EX_glc__D_e',
                     'EX_gln__L_e', 'EX_glu__L_e', 'EX_h_e', 'EX_h2o_e', 'EX_lac__D_e', 'EX_mal__L_e', 'EX_nh4_e',
                     'EX_o2_e',
                     'EX_pi_e', 'EX_pyr_e', 'EX_succ_e', 'FBA', 'FBP', 'FORt2', 'FORt', 'FRD7', 'FRUpts2', 'FUM',
                     'FUMt2_2',
                     'G6PDH2r', 'GAPD', 'GLCpts', 'GLNS', 'GLNabc', 'GLUDy', 'GLUN', 'GLUSy', 'GLUt2r', 'GND', 'H2Ot',
                     'ICDHyr', 'ICL', 'LDH_D', 'MALS', 'MALt2_2', 'MDH', 'ME1', 'ME2', 'NADH16', 'NADTRHD', 'NH4t',
                     'O2t',
                     'PDH', 'PFK', 'PFL', 'PGI', 'PGK', 'PGL', 'PGM', 'PIt2r', 'PPC', 'PPCK', 'PPS', 'PTAr', 'PYK',
                     'PYRt2',
                     'RPE', 'RPI', 'SUCCt2_2', 'SUCCt3', 'SUCDi', 'SUCOAS', 'TALA', 'THD2', 'TKT1', 'TKT2', 'TPI']

    # Dictionary to save the original bounds:
    original_bounds = {}

    # Cycle through the list and save the bounds in the dictionary:
    for reaction_id in reaction_list:
        reaction = model.reactions.get_by_id(reaction_id)
        original_bounds[reaction_id] = {
            'lower_bound': reaction.lower_bound,
            'upper_bound': reaction.upper_bound
        }

    # Display the dictionary with the original bounds:
    print("Original Bounds:", original_bounds)

    print(
        f"Original Bounds for Glucose: Lower = {original_bounds['EX_glc__D_e']['lower_bound']}, Upper = {original_bounds['EX_glc__D_e']['upper_bound']}")
    print(
        f"Original Bounds for Oxygen: Lower = {original_bounds['EX_o2_e']['lower_bound']}, Upper = {original_bounds['EX_o2_e']['upper_bound']}")
    print(
        f"Original Bounds for Biomass: Lower = {original_bounds['BIOMASS_Ecoli_core_w_GAM']['lower_bound']}, Upper = {original_bounds['BIOMASS_Ecoli_core_w_GAM']['upper_bound']}")

    # Listing of all exchange reactions:
    exchange_reactions = [rxn.id for rxn in model.exchanges]
    print(f"Exchange Reactions: {exchange_reactions}")

    return original_bounds

    # FBA for different glucose values:


def fba_different_glucose_values(model: Model):
    model = model.copy()
    # DataFrame for storing the results:
    results_df = pd.DataFrame(columns=['Glucose_Lower_Bound', 'Objective_Value', 'Status'])

    # List of glucose values:
    glucose_values = np.arange(-15, 0.1, 0.1)

    # Initialize the list of DataFrames to be merged:
    df_list = []

    # Set optimization goal back to biomass production:
    model.objective = 'BIOMASS_Ecoli_core_w_GAM'

    # Cycle through the list of glucose values:
    for glucose_value in glucose_values:
        # Setting the lower bounds for glucose:
        model.reactions.EX_glc__D_e.lower_bound = glucose_value

        # Performing the FBA:
        solution = model.optimize()

        # Create a temporary DataFrame for the current results:
        temp_df = pd.DataFrame({
            'Glucose_Lower_Bound': [glucose_value],
            'Objective_Value': [solution.objective_value],
            'Status': [solution.status]
        })

        df_list.append(temp_df)

    # Merge the DataFrames:
    results_df = pd.concat(df_list, ignore_index=True)

    return results_df


def standardized_fba_different_glucose_values_df(model: Model) -> pd.DataFrame:
    model = model.copy()

    glucose_values = np.arange(-15, 0.1, 0.1)
    df_list = []

    model.objective = 'BIOMASS_Ecoli_core_w_GAM'
    biomass_reaction = model.reactions.get_by_id('BIOMASS_Ecoli_core_w_GAM')
    total_carbon = abs(calculate_total_carbon_in_biomass(biomass_reaction))

    for glucose_value in glucose_values:
        model.reactions.get_by_id('EX_glc__D_e').lower_bound = glucose_value
        solution = model.optimize()

        if solution.status == 'optimal':
            temp_df = pd.DataFrame({
                'Glucose_Lower_Bound': [glucose_value*6],
                'Standardized_Biomass': [solution.objective_value*total_carbon],
                'Verhältnis': [(solution.objective_value*total_carbon)/abs(glucose_value*6)],
            })
            df_list.append(temp_df)

    results_df = pd.concat(df_list, ignore_index=True)
    return results_df

def calculate_total_carbon_in_biomass(biomass_reaction):
    lookup = get_carbon_count_lookup_biomass()
    total_carbon = 0
    for metabolite, coeff in biomass_reaction.metabolites.items():
        carbon_count = lookup.get(metabolite.id, 0)
        total_carbon += carbon_count * coeff
    return total_carbon


def visualize_biomass_vs_glucose(fba_different_glucose_values_df) -> None:
    # Creating a plot:
    plt.figure(figsize=(10, 6))

    # Inserting the data into the plot:
    plt.scatter(fba_different_glucose_values_df['Glucose_Lower_Bound'],
                fba_different_glucose_values_df['Objective_Value'], label='Biomass production')

    # Adding axis titles and a diagram title:
    plt.xlabel('Glucose Lower Bound')
    plt.ylabel('Biomass production')
    plt.title('Biomass production as a function of the glucose concentration')

    # Adding a legend:
    plt.legend()
    print(f"Saving biomass_vs_glucose.png")
    plt.savefig("depictions/biomass_vs_glucose.png")
    plt.show()

    # ------

    # FBA for different oxygen values:


def visualize_standardized_fba_different_glucose_values_df(fba_different_glucose_values_df) -> None:
    # Creating a plot:
    plt.figure(figsize=(10, 6))

    # Inserting the data into the plot:
    plt.scatter(fba_different_glucose_values_df['Glucose_Lower_Bound'],
                fba_different_glucose_values_df['Standardized_Biomass'], label='Biomass production')

    # Adding axis titles and a diagram title:
    plt.xlabel('Glucose Lower Bound [Bound*C]')
    plt.ylabel('Biomass production [Biomass*C]')
    plt.title('Standardized Biomass production as a function of the glucose concentration')

    # Adding a legend:
    plt.legend()
    print(f"Saving standardized_biomass_vs_glucose.png")
    plt.savefig("depictions/standardized_biomass_vs_glucose.png")
    plt.show()

    # ------

    # FBA for different oxygen values:


def fba_different_oxygen_values(model: Model, original_bounds: dict):
    # List of oxygen values:
    oxygen_values = np.arange(-100, 1, 0.1)

    # Initialize the list of DataFrames to be merged:
    df_list = []

    # Goal of optimization set back to biomass production:
    model.objective = 'BIOMASS_Ecoli_core_w_GAM'

    # Scroll through the list of oxygen values:
    for oxygen_value in oxygen_values:
        # Setting the lower bound for oxygen:
        model.reactions.EX_o2_e.lower_bound = oxygen_value

        # Performing the FBA:
        solution = model.optimize()

        # Create a temporary DataFrame for the current results:
        temp_df = pd.DataFrame({
            'Oxygen_Lower_Bound': [oxygen_value],
            'Objective_Value': [solution.objective_value],
            'Status': [solution.status]
        })

        df_list.append(temp_df)

    # Merge the DataFrames:
    results_df = pd.concat(df_list, ignore_index=True)

    # Reset the bounds for oxygen to the original values:
    model.reactions.EX_o2_e.lower_bound = original_bounds['EX_o2_e']['lower_bound']
    model.reactions.EX_o2_e.upper_bound = original_bounds['EX_o2_e']['upper_bound']

    # Save the DataFrame as a CSV file:
    # results_df.to_csv("fba_results_for_various_oxygen_values.csv", index=False)
    return results_df


def standardized_fba_different_oxygen_values_df(model: Model) -> pd.DataFrame:
    oxygen_values = np.arange(-100, 1, 0.1)
    df_list = []

    biomass_reaction = model.reactions.get_by_id('BIOMASS_Ecoli_core_w_GAM')
    total_carbon = abs(calculate_total_carbon_in_biomass(biomass_reaction))

    # Loop through the list of oxygen values:
    for oxygen_value in oxygen_values:
        model.reactions.get_by_id('EX_o2_e').lower_bound =oxygen_value
        solution = model.optimize()

        if solution.status == 'optimal':
            temp_df = pd.DataFrame({
                'Oxygen_Lower_Bound': [oxygen_value],
                'Standardized_Biomass': [solution.objective_value*total_carbon],
                'Ratio': [(solution.objective_value * total_carbon) / abs(oxygen_value * 6)],
            })
            df_list.append(temp_df)

    results_df = pd.concat(df_list, ignore_index=True)
    return results_df


def visualize_biomass_vs_oxygen(fba_different_oxygen_values_df) -> None:
    # Creating a plot:
    plt.figure(figsize=(10, 6))

    # Inserting the data into the plot:
    plt.scatter(fba_different_oxygen_values_df['Oxygen_Lower_Bound'], fba_different_oxygen_values_df['Objective_Value'],
                label='Biomass production')

    # Adding axis titles and a diagram title:
    plt.xlabel('Oxygen Lower Bound')
    plt.ylabel('Biomass production')
    plt.title('Biomass production as a function of the oxygen concentration')

    # Adding a legend:
    plt.legend()

    print(f"Saving biomass_vs_oxygen.png")
    plt.savefig("depictions/biomass_vs_oxygen.png")
    plt.show()


def visualize_standardized_fba_different_oxygen_values_df(fba_different_oxygen_values_df) -> None:
    # Creating a plot:
    plt.figure(figsize=(10, 6))

    # Inserting the data into the plot:
    plt.scatter(fba_different_oxygen_values_df['Oxygen_Lower_Bound'],
                fba_different_oxygen_values_df['Standardized_Biomass'], label='Biomass production')

    # Adding axis titles and a diagram title:
    plt.xlabel('Oxygen Lower Bound')
    plt.ylabel('Biomass production [Biomass*C]')
    plt.title('Standardized Biomass production as a function of the oxygen concentration')

    # Adding a legend:
    plt.legend()
    print(f"Saving standardized_biomass_vs_oxygen.png")
    plt.savefig("depictions/standardized_biomass_vs_oxygen.png")
    plt.show()


def update_exchange_fluxes(model: Model) -> Model:
    # All exchange fluxes = 0 whose metabolite in the ecm tool = 0:

    # List of metabolites that were used in the ecm tool:
    ecm_metabolites = ["M_ac_e", "M_acald_e", "M_akg_e", "M_co2_e", "M_etoh_e", "M_for_e", "M_fru_e", "M_fum_e",
                       "M_glc__D_e", "M_gln__L_e", "M_glu__L_e", "M_h2o_e", "M_lac__D_e", "M_mal__L_e", "M_nh4_e",
                       "M_o2_e", "M_pi_e", "M_pyr_e", "M_succ_e", "M_h_e"]

    # Loop that goes through all the exchange reactions in the model:
    for reaction in model.exchanges:
        # Verify that the metabolite was used in the reaction in the ecm tool:
        metabolite_in_ecm = any(met.id in ecm_metabolites for met in reaction.metabolites.keys())

        # Verify that the metabolite is in the model (to rule out new reactions):
        metabolite_in_model = any(met.id in model.metabolites for met in reaction.metabolites.keys())

        # If the metabolite was used in the model but not in the ecm tool, the lower and upper bounds are set to 0:
        if metabolite_in_model and not metabolite_in_ecm:
            reaction.lower_bound = 0
            reaction.upper_bound = 0

    return model
    # ----


def restrict_glucose_flow(model: Model) -> pd.DataFrame:
    print('Set all exchange fluxes to 0 where the metabolite in the ecm tool is 0...')
    model = update_exchange_fluxes(model.copy())

    # Restrict the flow for the C source (glucose):
    model.reactions.EX_glc__D_e.lower_bound = -1
    model.reactions.EX_glc__D_e.upper_bound = 0

    # Leave the other metabolites open:
    for reaction in [model.reactions.EX_o2_e, model.reactions.EX_h2o_e, model.reactions.EX_co2_e,
                     model.reactions.EX_h_e]:
        if reaction.id != "EX_glc__D_e":
            reaction.lower_bound = -1000
            reaction.upper_bound = 1000

    model.objective = 'ATPM'
    model.reactions.ATPM.lower_bound = 0
    model.reactions.ATPM.upper_bound = 1000
    solution = model.optimize()

    print(f"Optimal flow through the ATP-consuming reaction: {solution.objective_value}")

    # Save result in DataFrame:
    atp_results_df = pd.DataFrame({
        'Reaction': ['ATPM'],
        'Optimal_Flow': [solution.objective_value]
    })
    return atp_results_df


def calculate_max_and_max_standardized_ATP_for_every_reaction(model: Model,
                                                              core_conversions_df: pd.DataFrame) -> pd.DataFrame:
    model = model.copy()

    # Remove all elements where objective !== 0
    core_conversions_df = core_conversions_df[core_conversions_df["objective"] == 0]

    # Create DataFrame for the results:
    results = []

    # Iterate through each row of the ecm tool table:
    for index, row in core_conversions_df.iterrows():
        # Set metabolites:
        carbon_count = 0
        for metabolite, value in row.items():
            carbon_count += carbon_count_for_metabolite(metabolite, value)

            reaction_name = metabolite_to_reaction_name(metabolite)
            if not hasattr(model.reactions, reaction_name):
                continue
            reaction_obj = getattr(model.reactions, reaction_name)

            if abs(value) < 1e-10:  # oder abs(value) <= 1e-10:
                # print(f"Value = 0 or very small | Value: {value} | Set bounds for {reaction_name}: to 0|0
                reaction_obj.lower_bound = 0
                reaction_obj.upper_bound = 0
            elif value > 0:
                # print(f"Value > 0 | Value: {value} | Set bounds for {reaction_name}: to 0|{value}
                reaction_obj.lower_bound = 0
                reaction_obj.upper_bound = value
            else:
                # print(f"Value < 0 | Value: {value} | Set bounds for {reaction_name}: to {value}|0
                reaction_obj.lower_bound = value
                reaction_obj.upper_bound = 0

        carbon_count = carbon_count / 2

        # Calculate maximum ATP:
        model.objective = 'ATPM'
        model.reactions.ATPM.lower_bound = 0
        model.reactions.ATPM.upper_bound = 1000
        solution = model.optimize()

        atp_per_c = solution.objective_value / carbon_count

        # Save results:
        results.append({
            'Row': index,
            'Optimal_ATP': solution.objective_value,
            'Fluxes': {rxn.id: rxn.flux for rxn in model.exchanges},
            'ATP_per_C': atp_per_c,
            'carbon_count': carbon_count,
        })
    # Save results as a CSV file:
    results_df = pd.DataFrame(results)
    return results_df


def calculate_max_and_max_standardized_biomass_for_every_reaction(model: Model,
                                                                  core_conversions_df: pd.DataFrame) -> pd.DataFrame:
    model = model.copy()
    # Remove all elements where objective !== 0
    core_conversions_df = core_conversions_df[core_conversions_df["objective"] == 0]

    # Create DataFrame for the results:
    results = []

    # Iterate through each row of the ecm tool table:
    for index, row in core_conversions_df.iterrows():
        # Set metabolites:
        carbon_count = 0
        for metabolite, value in row.items():
            carbon_count += carbon_count_for_metabolite(metabolite, value)

            reaction_name = metabolite_to_reaction_name(metabolite)
            if not hasattr(model.reactions, reaction_name):
                continue
            reaction_obj = getattr(model.reactions, reaction_name)

            if abs(value) < 1e-10:
                # print(f"Value = 0 or very small | Set bounds for {reaction_name}: to 0|0")
                reaction_obj.lower_bound = 0
                reaction_obj.upper_bound = 0
            elif value > 0:
                # print(f"Value > 0 | Set bounds for {reaction_name}: to 0|{value}")
                reaction_obj.lower_bound = 0
                reaction_obj.upper_bound = value
            else:
                # print(f"Value < 0 | Set bounds for {reaction_name}: to {value}|0")
                reaction_obj.lower_bound = value
                reaction_obj.upper_bound = 0

        carbon_count = carbon_count / 2

        # Calculate maximum ATP:
        model.objective = 'BIOMASS_Ecoli_core_w_GAM'
        solution = model.optimize()
        atp_per_c = solution.objective_value / carbon_count

        # Save results:
        results.append({
            'Row': index,
            'Optimal_Biomass': solution.objective_value,
            'Fluxes': {rxn.id: rxn.flux for rxn in model.exchanges},
            'Biomass_per_C': atp_per_c,
            'carbon_count': carbon_count,
        })

    # Save results as a CSV file:
    results_df = pd.DataFrame(results)
    return results_df


def visualize_standardized_max_ATP(max_ATP_df: pd.DataFrame) -> None:
    equations = create_reaction_equations_atp_biomass(
        max_ATP_df.sort_values('ATP_per_C', ascending=False, ignore_index=True),
        'Fluxes', [0, 1, 2, 3, 4, 5], True)

    legend_elements = [Line2D([0], [0], color='w', marker='o', markerfacecolor='w', markersize=10, label=sf) for sf in
                       equations]

    max_ATP_df = max_ATP_df.sort_values(by='ATP_per_C', ascending=False)

    plt.figure(figsize=(14, 7))
    barWidth = 0.3

    # Use of absolute values for Gibbs Energy:
    plt.bar(range(len(max_ATP_df)), max_ATP_df['ATP_per_C'], color='r', label='Standardized maximum ATP')

    # Adding labels for the bars:
    plt.xlabel('Reactions', fontweight='bold')
    plt.ylabel('ATP/C')
    plt.title('Standardized maximum ATP values for the catabolic reactions of the E. coli core model')
    plt.xticks(range(len(max_ATP_df)), range(len(max_ATP_df)))

    plt.legend(handles=legend_elements, loc="upper right")

    # Save diagram:
    print(f"Saving standardized_max_ATP.png")
    plt.savefig('depictions/standardized_max_ATP.png')
    plt.show()


def visualize_standardized_max_biomass(max_biomass_df: pd.DataFrame) -> None:
    equations = create_reaction_equations_atp_biomass(
        max_biomass_df.sort_values('Biomass_per_C', ascending=False, ignore_index=True),
        'Fluxes', [0, 1, 2, 73, 74, 75], True)

    legend_elements = [Line2D([0], [0], color='w', marker='o', markerfacecolor='w', markersize=10, label=sf) for sf in
                       equations]

    max_biomass_df = max_biomass_df.sort_values(by='Biomass_per_C', ascending=False)

    plt.figure(figsize=(14, 7))
    barWidth = 0.3

    # Use of absolute values for Gibbs Energy:
    plt.bar(range(len(max_biomass_df)), max_biomass_df['Biomass_per_C'], color='r',
            label='Standardized maximum biomass')

    # Adding labels for the bars:
    plt.xlabel('Reactions', fontweight='bold')
    plt.ylabel('Biomass/C')
    plt.title('Standardized maximum biomass values for the catabolic reactions of the E. coli core model')
    plt.xticks(range(len(max_biomass_df)), range(len(max_biomass_df)))

    plt.legend(handles=legend_elements, loc="lower left")

    # Save diagram:
    print(f"Saving standardized_max_biomass.png")
    plt.savefig('depictions/standardized_max_biomass.png')
    plt.show()
