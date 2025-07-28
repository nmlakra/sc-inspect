# sc-inspect: A Terminal UI for Single-Cell Data Exploration

**sc-inspect** is a fast, lightweight, and interactive Terminal User Interface (TUI) for exploring the metadata of single-cell datasets stored in `.h5ad` format.

*(A placeholder for a cool GIF or screenshot of your application in action!)*

## 💡 The Problem

Exploring single-cell RNA-seq data often starts with examining the `AnnData` object, particularly the observation metadata (`adata.obs`). This usually involves firing up a Jupyter notebook, writing boilerplate `pandas` and `matplotlib` code just to get a basic sense of the available annotations, their data types, and distributions. This process can be repetitive and slow.

`sc-inspect` was built to solve this. It provides a zero-boilerplate, in-terminal way to get a quick overview of your `.h5ad` file's metadata, letting you inspect and visualize key information instantly.

## ✨ Features

* **Interactive Metadata Table**: View your `adata.obs` dataframe in a rich, scrollable table.

* **Column Details**: Select a column to view its name, data type, and summary statistics.

* **Instant Plotting**:

  * **Histograms**: Automatically generates a terminal-based histogram for any selected numerical column.

  * **UMAP Plots**: Visualizes `adata.obsm['X_umap']` colored by the selected categorical variable.

* **Built-in Demo**: Comes with a pre-packaged `pbmc3k` dataset so you can try it out immediately.

* **Fast & Responsive**: Built with the modern [Textual](https://github.com/Textualize/textual) TUI framework for a smooth experience.

## 🛠️ How It's Built

`sc-inspect` is built entirely in Python using a suite of powerful libraries:

* [**Textual**](https://github.com/Textualize/textual): For the core TUI application framework, layout, and widgets.

* [**Rich**](https://github.com/Textualize/rich): For beautiful and informative rendering of tables, logs, and plots in the terminal.

* [**plotille**](https://github.com/tammoippen/plotille): For generating lightweight, text-based plots directly in the terminal.

* [**Scanpy**](https://scanpy.readthedocs.io/en/latest/): For reading `.h5ad` files and handling single-cell data structures.

## 🚀 Getting Started

### Prerequisites

* Python 3.10+

* [uv](https://github.com/astral-sh/uv) (recommended) or `pip`

### Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/nmlakra/sc-inspect.git
   cd sc-inspect
   ```

2. Create and activate a virtual environment:

   Using uv:
   ```bash
   uv venv
   source .venv/bin/activate
   ```

   Using pip:
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   ```

3. Install the package:

   Using uv:
   ```bash
   uv pip install .
   ```

   Using pip:
   ```bash
   pip install .
   ```

### Usage

The `sc-inspect` command will be available automatically after installation.

**1. Run with the built-in demo file:**

Simply run the command without any arguments to explore the included `pbmc3k` dataset.

```bash
sc-inspect
```

**2. Run with your own file:**

Provide the path to your `.h5ad` file.

```bash
sc-inspect /path/to/your/data.h5ad
```

## ⌨️ How to Use

* **Navigate Table**: Use the **Up** and **Down** arrow keys to select different metadata columns in the left-hand table.

* **View Plots**: As you select a column, a relevant plot (histogram for numeric, UMAP for categorical) will automatically appear in the bottom-right panel.

* **Quit**: Press **`q`** at any time to exit the application.

## 🧑‍💻 For Developers

To set up a development environment, install the package in editable (`-e`) mode with the `dev` dependencies:

```bash
uv pip install -e .[dev]
```

This will install development tools like `black` and `isort`.


## 🔮 What's Next for sc-inspect

* \[ \] Support for viewing `adata.var` and exploring gene expression in the data.

* \[ \] Add more plot types (e.g., violin plots for numeric data).

* \[ \] Implement search/filtering for the metadata table.

* \[ \] Allow customization of UMAP plots (e.g., changing plot size, color palettes).

* \[ \] Package the project for easy distribution via `pip`.

