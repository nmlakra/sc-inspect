import scanpy as sc
import plotille
from rich.text import Text

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


def get_umap_key(adata: sc.AnnData) -> str:
    """To avoid cases with different casing of the key: 'X_umap', 'X_Umap', 'X_UMAP', etc."""
    umap_key = ""

    for key in adata.obsm_keys():
        if "umap" in key.lower():
            umap_key = key

    return umap_key


def plot_umap(adata: sc.AnnData, column) -> Text:
    umap_key = get_umap_key(adata)

    if not umap_key:
        return Text(
            "Error: cannot plot UMAP, adata doesn't contain the UMAP information."
        )

    if column not in adata.obs.columns:
        return Text("Error: adata invalid column selected.")

    fig = plotille.Figure()
    fig.background = 232
    fig._origin = False
    fig.width = 40
    fig.height = 20
    fig.color_mode = "byte"

    for i, val in enumerate(adata.obs[column].unique()):
        bdata = adata[adata.obs[column] == val]

        x = bdata.obsm["X_umap"][:, 0]
        y = bdata.obsm["X_umap"][:, 1]

        label_color = plotille._colors.rgb2byte(*(COLORS[i % len(COLORS)]))

        fig.scatter(x, y, label=val, lc=label_color)

    umap = remove_umap_spline(fig.show(legend=True))
    return Text.from_ansi(umap)


def plot_hist(adata: sc.AnnData, column) -> Text:
    data_series = adata.obs[column]
    hist = plotille.hist(
        data_series, bg=232, color_mode="byte", lc=231, width=25, bins=25
    )
    return Text.from_ansi(hist)


def remove_umap_spline(umap_fig: str) -> str:
    cleaned_umap = []
    lines = umap_fig.splitlines()
    offset = lines[1].find(" |") + 3
    for i, line in enumerate(lines):
        if line == "":
            cleaned_umap.extend(lines[i:])
            break
        elif i == 0:
            continue
        elif line.endswith("(X)"):
            continue
        elif line.startswith(" "):
            continue
        else:
            cleaned_umap.append(" " * offset + line[offset:])
    return "\n".join(cleaned_umap)
