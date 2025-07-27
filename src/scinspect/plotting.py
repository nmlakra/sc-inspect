import scanpy as sc
import plotille

COLORS = [
    (166, 206, 227),  # Sky
    (31, 120, 180),  # Blue
    (178, 223, 138),  # Lime
    (51, 160, 44),  # Green
    (251, 154, 153),  # Salmon
    (227, 26, 28),  # Red
    (253, 191, 111),  # Peach
    (255, 127, 0),  # Orange
    (202, 178, 214),  # Lilac
    (106, 61, 154),  # Purple
    (255, 255, 153),  # Cream
    (177, 89, 40),  # Brown
]


def plot_umap(adata: sc.AnnData, column) -> str:
    if "X_umap" not in adata.obsm:
        return "Error: cannot plot UMAP, adata doesn't contain the UMAP information."

    if column not in adata.obs.columns:
        return "Error: adata invalid column selected."

    fig = plotille.Figure()
    fig._origin = False
    fig.width = 30
    fig.height = 15
    fig.color_mode = "rgb"

    for i, val in enumerate(adata.obs[column].unique()):
        bdata = adata[adata.obs[column] == val]

        x = bdata.obsm["X_umap"][:, 0]
        y = bdata.obsm["X_umap"][:, 1]

        label_color = COLORS[i % len(COLORS)]

        fig.scatter(x, y, label=val, lc=label_color)

    return fig.show(legend=True)
