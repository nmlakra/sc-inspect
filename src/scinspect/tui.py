import scanpy as sc
from rich.align import Align
from rich.table import Table
from rich.text import Text
from textual.app import App, ComposeResult
from textual.containers import Horizontal
from textual.widgets import DataTable, Footer, Header, RichLog

from .plotting import plot_hist, plot_umap


def read_h5ad(file_path: str) -> sc.AnnData:
    adata = sc.read_h5ad(file_path, backed="r+")
    return adata


class ScInspectApp(App):
    def __init__(self, filepath: str) -> None:
        self.filepath = filepath
        self.adata = None
        super().__init__()

    CSS_PATH = "scinspect.tcss"

    def compose(self) -> ComposeResult:
        yield Header()
        yield Footer()
        with Horizontal(id="main-container"):
            yield DataTable(id="obs-table", cursor_type="row")

            with Horizontal(id="details-container"):
                yield RichLog(id="details-log")
                yield RichLog(id="plot-log")

    @staticmethod
    def translate_dtype(dtype: str) -> str:
        if not dtype:
            return "unkown"

        translation_map = {
            "int32": "numeric",
            "int64": "numeric",
            "float32": "numeric",
            "float64": "numeric",
            "category": "categorical",
        }

        if dtype in translation_map:
            return translation_map[dtype]
        return dtype

    def get_data_columns(self) -> list:
        data_list = [("Column Name", "Data Type")]
        if self.adata is None:
            raise Exception("could not load the data")
        # total_counts = len(self.adata.obs)
        for col in self.adata.obs.columns:
            # non_nan_counts = str(total_counts - sum(self.adata.obs[col].isna()))
            col_dtype = self.translate_dtype(str(self.adata.obs[col].dtype))
            data_list.append((str(col), col_dtype))

        return data_list

    def on_mount(self) -> None:
        self.title = "ScInspect"
        self.sub_title = self.filepath

        self.adata = read_h5ad(self.filepath)

        plot_log = self.query_one("#plot-log", RichLog)
        plot_log.border_title = "Plot"

        details_log = self.query_one("#details-log", RichLog)
        details_log.border_title = "Data Summary"

        if self.adata:
            rows = self.get_data_columns()
            table = self.query_one("#obs-table", DataTable)
            table.add_columns(*rows[0])
            for idx, row in enumerate(rows[1:], start=1):
                label = Text(str(idx), style="#B0FC38")
                table.add_row(*row, label=label, key=row[0])

    def on_data_table_row_selected(self, event: DataTable.RowSelected) -> None:
        if self.adata:
            details_log = self.query_one("#details-log", RichLog)
            plot_log = self.query_one("#plot-log", RichLog)

            selected_row = event.row_key.value
            if selected_row is None:
                raise Exception("Error: selected row doesn't have a vaild label")

            details_log.clear()
            plot_log.clear()
            data_series = self.adata.obs[selected_row]
            col_type = str(data_series.dtype)

            if "category" in col_type:
                details_log.write(self.get_categorical_summary_table(selected_row))
                details_log.write(f"  Unique Values: {data_series.nunique()}")
                details_log.write(f"  Total Counts: {len(data_series)}")

                plot_log.write(plot_umap(self.adata, selected_row))

            elif "int" in col_type or "float" in col_type:
                details_log.write(self.get_numeric_summary_table(selected_row))

                plot_log.write(Align.center(plot_hist(self.adata, selected_row)))

    def get_categorical_summary_table(self, categorical_col: str) -> Table:
        # TODO: add NA count and NA%
        summary_table = Table(title=categorical_col)

        if self.adata is None:
            raise Exception("Error: data is not loaded")
        if categorical_col not in self.adata.obs.columns:
            raise Exception(f"Error: {categorical_col} not found in the .h5ad file")

        summary_table.add_column("Value")
        summary_table.add_column("Count", justify="right")

        data_series = self.adata.obs[categorical_col]
        for row in data_series.value_counts().items():
            summary_table.add_row(*(str(data) for data in row))

        return summary_table

    def get_numeric_summary_table(self, numeric_col: str) -> Table:
        summary_table = Table(title=numeric_col)

        if self.adata is None:
            raise Exception("Error: data is not loaded")
        if numeric_col not in self.adata.obs.columns:
            raise Exception(f"Error: {numeric_col} not found in the .h5ad file")

        summary_table.add_column("Stat")
        summary_table.add_column("Value", justify="right")

        data_series = self.adata.obs[numeric_col]
        data_count = data_series.count()
        data_mean = data_series.mean()
        data_std = data_series.std()
        data_min = data_series.min()
        data_lower_quartile = data_series.quantile(0.25)
        data_median = data_series.median()
        data_upper_quartile = data_series.quantile(0.75)
        data_max = data_series.max()
        data_na_count = data_series.isna().sum()
        data_na_percent = (data_na_count / data_count) * 100

        summary_table.add_row("Count", str(data_count))
        summary_table.add_row("Mean", str(data_mean))
        summary_table.add_row("SD", str(data_std))
        summary_table.add_row("Min", str(data_min))
        summary_table.add_row("25%", str(data_lower_quartile))
        summary_table.add_row("50%", str(data_median))
        summary_table.add_row("75%", str(data_upper_quartile))
        summary_table.add_row("Max", str(data_max))
        summary_table.add_row("Null Count", str(data_na_count))
        summary_table.add_row("Null%", str(data_na_percent))

        return summary_table
