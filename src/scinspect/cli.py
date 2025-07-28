import argparse
from importlib import resources

from .tui import ScInspectApp


def main() -> None:
    """Entry point for the sc-inspect command."""
    parser = argparse.ArgumentParser(
        description="A TUI to inspect single-cell data.",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "file",
        nargs="?",
        help="Path to the .h5ad file. If not provided, the provided demo file is used.",
    )
    args = parser.parse_args()

    if args.file:
        filepath = args.file
    else:
        with resources.path("scinspect.examples", "pbmc3k_processed.h5ad") as path:
            filepath = str(path)
            print(f"No file path was provided. Using the demo file: {filepath}")

    app = ScInspectApp(filepath=filepath)
    app.run()


if __name__ == "__main__":
    main()
