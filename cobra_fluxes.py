from multiprocessing import freeze_support

import matplotlib.pyplot as plt
import pandas as pd
from cobra import Reaction, Model
from matplotlib.lines import Line2D
from scripts.helper import get_carbon_count_lookup


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


def fba_different_glucose_values(model: Model, original_bounds: dict):
    # DataFrame for storing the results:
    results_df = pd.DataFrame(columns=['Glucose_Lower_Bound', 'Objective_Value', 'Status'])

    # List of glucose values:
    glucose_values = [-10, -8, -6, -5, -4, -2, -1, 0]

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

        # Verify that the solution is valid:
        if solution.status == 'infeasible':
            print(f"Warning: solution is not allowed for glucose value {glucose_value}.")
            continue

        # Create a temporary DataFrame for the current results:
        temp_df = pd.DataFrame({
            'Glucose_Lower_Bound': [glucose_value],
            'Objective_Value': [solution.objective_value],
            'Status': [solution.status]
        })

        df_list.append(temp_df)

    # Merge the DataFrames:
    results_df = pd.concat(df_list, ignore_index=True)

    # Reset the bounds for glucose to the original values:
    model.reactions.EX_glc__D_e.lower_bound = original_bounds['EX_glc__D_e']['lower_bound']
    model.reactions.EX_glc__D_e.upper_bound = original_bounds['EX_glc__D_e']['upper_bound']

    return results_df


def visualize_biomass_vs_glucose(fba_different_glucose_values_df) -> None:
    # Creating a plot:
    plt.figure(figsize=(10, 6))

    # Inserting the data into the plot:
    plt.scatter(fba_different_glucose_values_df['Glucose_Lower_Bound'],
                fba_different_glucose_values_df['Objective_Value'], label='Biomasseproduktion')

    # Adding axis titles and a diagram title:
    plt.xlabel('Glucose Lower Bound')
    plt.ylabel('Objective Value (Biomass production)')
    plt.title('Biomass production as a function of glucose concentration')

    # Adding a legend:
    plt.legend()
    print(f"Saving biomass_vs_glucose.png")
    plt.savefig("depictions/biomass_vs_glucose.png")
    plt.show()

    # ------

    # FBA for different oxygen values:


def fba_different_oxygen_values(model: Model, original_bounds: dict):
    # DataFrame for storing the results:
    results_df = pd.DataFrame(columns=['EX_o2_e_Lower_Bound', 'Objective_Value', 'Status'])

    # List of oxygen values:
    oxygen_values = [-30, -25, -20, -15, -10, -5, 0]

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

        # Verify that the solution is valid:
        if solution.status == 'infeasible':
            print(f"Warning: Solution is not allowed for oxygen value {oxygen_value}.")
            continue

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


def visualize_biomass_vs_oxygen(fba_different_oxygen_values_df) -> None:
    # Creating a plot:
    plt.figure(figsize=(10, 6))

    # Inserting the data into the plot:
    plt.scatter(fba_different_oxygen_values_df['Oxygen_Lower_Bound'], fba_different_oxygen_values_df['Objective_Value'],
                label='Biomass production')

    # Adding axis titles and a diagram title:
    plt.xlabel('Oxygen Lower Bound')
    plt.ylabel('Objective Value (Biomass production)')
    plt.title('Biomass production as a function of oxygen concentration')

    # Adding a legend:
    plt.legend()

    print(f"Saving biomass_vs_oxygen.png")
    plt.savefig("depictions/biomass_vs_oxygen.png")
    plt.show()

    # ------


def add_ATP_hydrolysis_reaction(model: Model) -> Model:
    # Create a new reaction:
    new_reaction = Reaction('ATP_hydrolysis')
    new_reaction.name = 'ATP Hydrolysis'

    # Defining the reaction:
    new_reaction.add_metabolites({
        model.metabolites.atp_c: -1,
        model.metabolites.h2o_c: -1,
        model.metabolites.adp_c: 1,
        model.metabolites.pi_c: 1  # or model.metabolites.pi_e, if the reaction takes place in the extracellular space
    })

    # Adding the reaction to the model:
    model.add_reactions([new_reaction])
    return model
    # ----


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
    # Restrict the flow for the C source (glucose):
    model.reactions.EX_glc__D_e.lower_bound = -1
    model.reactions.EX_glc__D_e.upper_bound = 0

    # Leave the other metabolites open:
    for reaction in [model.reactions.EX_o2_e, model.reactions.EX_h2o_e, model.reactions.EX_co2_e,
                     model.reactions.EX_h_e]:
        if reaction.id != "EX_glc__D_e":
            print(reaction.id)
            reaction.lower_bound = -1000
            reaction.upper_bound = 1000

    model.objective = 'ATP_hydrolysis'
    solution = model.optimize()
    print(f"Optimal flow through the ATP-consuming reaction: {solution.objective_value}")

    # Save result in DataFrame:
    atp_results_df = pd.DataFrame({
        'Reaction': ['ATP_hydrolysis'],
        'Optimal_Flow': [solution.objective_value]
    })
    return atp_results_df


def temp(model, original_bounds):
    # Reset the bounds for the metabolites to the original values: (?)
    for reaction_id, bounds in original_bounds.items():
        model.reactions.get_by_id(reaction_id).lower_bound = bounds['lower_bound']
        model.reactions.get_by_id(reaction_id).upper_bound = bounds['upper_bound']

    # ------

    # Lookup:
    carbon_count_lookup = get_carbon_count_lookup()

    # ecm Tool read in table:
    ecm_data = pd.read_csv("csvs/core_conversions.csv")
    ecm_data = ecm_data[ecm_data["objective"] == 0]

    # Create DataFrame for the results:
    results = []

    print(model.reactions.get_by_id("EX_ac_e"))

    # Function to convert metabolite name from table to reaction name in model:
    def metabolite_to_reaction_name(metabolite_name):
        if metabolite_name == "objective":
            return "Biomass_Ecoli_core"
        else:
            return "EX_" + metabolite_name.split('_')[1] + "_e"

    def carbon_count_for_metabolite(name, val):
        multiplier = carbon_count_lookup.get(name, 0)
        return abs(val) * multiplier

    print('-----------------------------------------------')
    atp_per_c = []
    # Iterate through each row of the ecm tool table:
    for index, row in ecm_data.iterrows():
        # Set metabolites:
        carbon_count = 0
        for metabolite, value in row.items():
            carbon_count += carbon_count_for_metabolite(metabolite, value)

            reaction_name = metabolite_to_reaction_name(metabolite)
            if not hasattr(model.reactions, reaction_name):
                continue
            reaction_obj = getattr(model.reactions, reaction_name)

            if value == 0:
                # print(f"Value = 0 | Set bounds for {reaction_name}: to 0|0")
                reaction_obj.lower_bound = 0
                reaction_obj.upper_bound = 0
            elif value > 0:
                # print(f"Value > 0 | Set bounds for {reaction_name}: to 0|1000")
                reaction_obj.lower_bound = 0
                reaction_obj.upper_bound = 1000
            else:
                # print(f"Value < 0 | Set bounds for {reaction_name}: to {value}|0")
                reaction_obj.lower_bound = value
                reaction_obj.upper_bound = 0

        carbon_count = carbon_count / 2

        # Calculate maximum ATP:
        model.objective = 'ATP_hydrolysis'
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
    # exit()

    # Save results as a CSV file:
    results_df = pd.DataFrame(results)
    results_df.to_csv("max_ATP_for_every_reaction.csv", index=False)

    # Plotting:
    summenformeln = ["Reaction 0: 10 C6H12O6 + 0.0 O2 + 1.2 PO₄³⁻ -> 0.0 H2O + 20 C3H6O3 + 0.0 C4H6O4",
                     "Reaction 1: 1.4 CH2O2 + 10 C6H12O6 -> 20 C3H6O3",
                     "Reaction 2:1.4 CH2O2 + 10 C6H12O6 -> 20 C3H6O3",
                     "Reaction 3: 10 C6H12O6 + 0.1 O2 -> 0.2 H2O + 19.8 C3H6O3 + 0.2 C4H6O4",
                     "Reaction 4: 1.4 CH2O2 + 10 C6H12O6 -> 20 C3H6O3",
                     "Reaction 5: 10 C6H12O6 + 0.3 O2 -> 0.2 CO2 + 0.3 H2O + 19.9 C3H6O3 + 2.5 PO₄³⁻"]

    legend_elements = [Line2D([0], [0], color='w', marker='o', markerfacecolor='w', markersize=10, label=sf) for sf in
                       summenformeln]

    max_ATP_df = pd.read_csv("max_ATP_for_every_reaction.csv")

    max_ATP_df = max_ATP_df.sort_values(by='ATP_per_C', ascending=False)

    plt.figure(figsize=(14, 7))
    barWidth = 0.3

    # Use of absolute values for Gibbs Energy:
    plt.bar(range(len(max_ATP_df)), max_ATP_df['ATP_per_C'].abs(), color='r', label='Standardized maximum ATP')

    # Adding labels for the bars:
    plt.xlabel('Reactions', fontweight='bold')
    plt.ylabel('ATP/C')
    plt.title('Standardized maximum ATP values for the catabolic reactions of the E. coli core model')
    plt.xticks(range(len(max_ATP_df)), range(len(max_ATP_df)))

    plt.legend(handles=legend_elements, loc="upper right")

    # Save diagram:
    plt.savefig('Standardized_max_ATP.png')
    plt.show()
