from enum import Enum
from pymatgen.core.periodic_table import Element
import numpy as np

import plotly.graph_objects as go


class VisualizationSetting(Enum):
    NONE = 0
    LAST = 1
    ALL = 2
    ALL_DETAILED = 3


def element_color(atomic_number: int) -> str:
    # Dictionary mapping chemical symbols to colors
    color_map = {
        1: "#F0F8FF",  # Hydrogen
        2: "#D3D3D3",  # Helium
        3: "#B22222",  # Lithium
        4: "#9ACD32",  # Beryllium
        5: "#FFD700",  # Boron
        6: "#000000",  # Carbon
        7: "#8F00FF",  # Nitrogen
        8: "#FF0000",  # Oxygen
        9: "#DAA520",  # Fluorine
        10: "#EE82EE",  # Neon
        11: "#0000FF",  # Sodium
        12: "#228B22",  # Magnesium
        13: "#C0C0C0",  # Aluminum
        14: "#A52A2A",  # Silicon
        15: "#FFA500",  # Phosphorus
        16: "#FFFF00",  # Sulfur
        17: "#00FF00",  # Chlorine
        18: "#FF00FF",  # Argon
        19: "#F0E68C",  # Potassium
        20: "#E6E6FA",  # Calcium
    }

    # Return the color for the given chemical symbol, or a default color if not found
    return color_map.get(atomic_number, "#808080")  # Default color is gray


def visualize_and_save_crystal(
    atomic_numbers: np.ndarray,
    frac_x: np.ndarray,
    name: str,
):
    fig = plot_crystal(atomic_numbers, frac_x)
    # Save the plot as a PNG file
    fig.write_image(name + ".png")
    return fig


def plot_crystal(
    atomic_numbers: np.ndarray,
    coords: np.ndarray,
) -> go.Figure:
    # we use min(z, 118) since the mask state needs a visualization
    element_symbols = [Element.from_Z(min(z, 118)).symbol for z in atomic_numbers]
    atom_idx = np.arange(len(atomic_numbers)).tolist()
    pos_arr = coords
    # Create a Plotly figure
    fig = go.Figure()

    # fig.add_trace(go.Scatter3d(
    #     x=list(map(lambda x : x[0], pos_arr)),
    #     y=list(map(lambda x : x[1], pos_arr)),
    #     z=list(map(lambda x : x[2], pos_arr)),
    #     mode='markers',
    #     marker=dict(size=5,
    #         color=list(map(element_color, atomic_numbers)),  # Set the color based on the atom type
    #     ),
    #     text=element_symbols,
    #     name=atom_idx,
    #     hoverinfo='text+x+y+z+name',
    # ))

    # Add traces for each atom in the structure
    for i in range(len(pos_arr)):
        coord = pos_arr[i]
        atomic_num = atomic_numbers[i]
        fig.add_trace(
            go.Scatter3d(
                x=[coord[0]],
                y=[coord[1]],
                z=[coord[2]],
                mode="markers",
                marker=dict(
                    size=5,
                    color=element_color(
                        atomic_num
                    ),  # Set the color based on the atom type
                ),
                text=element_symbols[i],
                name=atom_idx[i],
            )
        )

    # Set the layout for the 3D plot
    fig.update_layout(
        title="Crystal Structure",
        scene=dict(xaxis_title="X", yaxis_title="Y", zaxis_title="Z"),
        margin=dict(l=0, r=0, b=0, t=0),
    )
    return fig