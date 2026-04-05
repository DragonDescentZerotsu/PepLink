from pathlib import Path


def main() -> None:
    notebook = Path(__file__).with_suffix(".ipynb")
    print(f"See {notebook.name} for the full set of runnable PepLink examples.")


if __name__ == "__main__":
    main()
