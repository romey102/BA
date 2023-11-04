import numpy as np
import pandas as pd
from equilibrator_api import ComponentContribution, Q_


def calculate_enthalpies_of_formation(core_conversions_df: pd.DataFrame) -> pd.DataFrame:
    # Extract catabolic reactions
    core_conversions_df = core_conversions_df[core_conversions_df["objective"] == 0]

    core_conversions_df = core_conversions_df.reset_index(drop=True)

    # there is no DGf for protons and objective (biomass)
    core_conversions_df = core_conversions_df.drop('M_h_e', axis=1)
    core_conversions_df = core_conversions_df.drop('objective', axis=1)

    # normalize
    core_conversions_df = core_conversions_df.div(core_conversions_df.abs().sum(axis=1), axis=0)

    # clean columns
    columns = list(core_conversions_df.columns)
    metabolite_names = [word.replace("M_", "").replace("_e", "") for word in columns]

    cc = ComponentContribution()
    # if you want to change the aqueous environment
    cc.p_h = Q_(7.4)
    cc.p_mg = Q_(3.0)
    cc.ionic_strength = Q_("0.25M")
    cc.temperature = Q_("298.15K")

    # obtain list of compound objects
    compound_list = [cc.get_compound(f"bigg.metabolite:{cname}") for cname in metabolite_names]

    # apply enthalpies_of_formation
    standard_dgf_mu, sigmas_fin, sigmas_inf = zip(*map(cc.standard_dg_formation, compound_list))
    standard_dgf_mu = np.array(standard_dgf_mu)

    core_conversions_df = pd.DataFrame()
    for name, value in zip(columns, standard_dgf_mu):
        core_conversions_df[name] = [value]

    return core_conversions_df
