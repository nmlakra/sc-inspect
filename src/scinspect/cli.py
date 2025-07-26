import argparse
from .tui import ScInspectApp


def main() -> None:
    parser = argparse.ArgumentParser(description="A TUI to inspect single-cell data.")
    parser.add_argument("file", help="Path to the .h5ad file")
    args = parser.parse_args()

    app = ScInspectApp(filepath=args.file)
    app.run()


if __name__ == "__main__":
    main()
