from textual.app import App, ComposeResult
from textual.widgets import Markdown, RichLog, Header
import scanpy as sc


def read_h5ad(file_path: str) -> sc.AnnData:
    adata = sc.read_h5ad(file_path)
    return adata


class ScInspectApp(App):
    def __init__(self, filepath: str) -> None:
        self.filepath = filepath
        self.adata = None
        super().__init__()

    def compose(self) -> ComposeResult:
        yield Header()
        yield RichLog()

    def on_mount(self) -> None:
        self.title = "ScInspect"
        self.sub_title = self.filepath

        self.adata = read_h5ad(self.filepath)

    def on_ready(self) -> None:
        text_log = self.query_one(RichLog)

        adata_info = self.adata.__repr__()
        text_log.write(adata_info)
