import os

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx

from scripts.helper import metabolite_name_lookup

DEPICTIONS_PATH = 'depictions/flux_network_reactions/'


def visualize_flux_network_reactions(flux_network_diagram_df: pd.DataFrame) -> None:
    os.makedirs(DEPICTIONS_PATH, exist_ok=True)

    # Catabolic reactions only:
    flux_network_diagram_df = flux_network_diagram_df[flux_network_diagram_df["objective"] == 0]

    # Lookup dictionary of metabolite names:
    lookup_dict = metabolite_name_lookup()

    # Colors:
    reactant_color = 'red'
    product_color = 'green'

    # Display each reaction in a separate diagram:
    for idx, row in flux_network_diagram_df.iterrows():
        plt.figure(figsize=(12, 12))
        plt.title('Flowchart for Individual Catabolic Reactions in the E. Coli Core Network')  # Headline
        G = nx.DiGraph()
        edge_colors = []
        reaction = f"Reaction_{idx + 1}"

        for metabolite, value in row.items():
            if value < 0:
                G.add_edge(metabolite, reaction, weight=-value)
                edge_colors.append(reactant_color)
            elif value > 0:
                G.add_edge(reaction, metabolite, weight=value)
                edge_colors.append(product_color)

        # Layout and drawing:
        pos = nx.spring_layout(G, seed=42)
        nx.draw(G, pos,
                edge_color=edge_colors,
                node_shape='o',  # Standard shape for all nodes
                node_size=700,
                font_size=18)

        # Legend for reactants and products:
        plt.plot([0], [0], color=reactant_color, label='Reactant')
        plt.plot([0], [0], color=product_color, label='Product')
        plt.legend(loc='upper right')

        # Add metabolite names:
        labels = {n: lookup_dict.get(n, n) for n in G.nodes()}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=10, font_color='black')

        # Edge Weight Labels:
        labels = {k: round(v, 3) for k, v in nx.get_edge_attributes(G, 'weight').items()}
        nx.draw_networkx_edge_labels(G, pos, edge_labels=labels, label_pos=0.7, font_size=10)

        # Add heading:
        plt.title('Flow Diagram for Individual Catabolic Reactions in E. Coli Core Network')

        # Save diagrams:
        file_name = f'flux_network_reaction_{idx + 1}.png'
        print(f"Saving {file_name}")
        plt.savefig(DEPICTIONS_PATH + file_name)
        plt.close()
