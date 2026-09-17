"""Extract GEO series-matrix expression tables for Galaxy input.

Features
--------
- Searches only directories whose names begin with GSE
- Searches recursively inside each GSE directory
- Reads *_series_matrix.txt files
- Reads *_series_matrix.txt.gz files directly
- Produces Galaxy-ready intensity CSV files
- Can also be imported and called from another Python pipeline

Default output:
    data/galaxy_file_input/
"""

import argparse
import csv
import gzip
from pathlib import Path


DEFAULT_OUTPUT_DIR = Path("data/galaxy_file_input")


def make_output_name(input_file: Path) -> str:
    """Create a consistent intensity CSV filename."""

    name = input_file.name

    if name.lower().endswith(".gz"):
        name = name[:-3]

    if name.lower().endswith("_series_matrix.txt"):
        name = name[:-len("_series_matrix.txt")]
    elif name.lower().endswith(".txt"):
        name = name[:-4]

    name = name.replace("-", "_")

    return f"{name}_intensities.csv"


def extract_matrix(lines, output_file: Path, source_name: str) -> Path:
    """Extract the GEO matrix between begin/end markers."""

    inside_matrix = False
    found_end = False
    rows_written = 0
    column_count = 0

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with output_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as dst:

        writer = csv.writer(dst)

        for line in lines:
            marker = line.strip()

            if marker == "!series_matrix_table_begin":
                inside_matrix = True
                continue

            if marker == "!series_matrix_table_end":
                found_end = True
                break

            if not inside_matrix:
                continue

            row = next(
                csv.reader(
                    [line],
                    delimiter="\t",
                    quotechar='"',
                )
            )

            if rows_written == 0:
                if row and row[0] == "ID_REF":
                    row[0] = "probe_id"

                column_count = len(row)

            writer.writerow(row)
            rows_written += 1

    if not inside_matrix:
        raise RuntimeError(
            f"!series_matrix_table_begin not found in {source_name}"
        )

    if not found_end:
        raise RuntimeError(
            f"!series_matrix_table_end not found in {source_name}"
        )

    if rows_written == 0:
        raise RuntimeError(
            f"No expression matrix found in {source_name}"
        )

    print()
    print("Matrix extraction complete")
    print(f"Source:         {source_name}")
    print(f"Output:         {output_file}")
    print(f"Probe rows:     {rows_written - 1}")
    print(f"Sample columns: {max(column_count - 1, 0)}")

    return output_file


def process_matrix_file(input_file: Path, output_dir: Path) -> Path:
    """Process either TXT or TXT.GZ GEO series matrix."""

    output_file = output_dir / make_output_name(input_file)

    if input_file.name.lower().endswith(".gz"):

        with gzip.open(
            input_file,
            "rt",
            encoding="utf-8",
            errors="replace",
            newline="",
        ) as src:

            return extract_matrix(
                src,
                output_file,
                str(input_file),
            )

    with input_file.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as src:

        return extract_matrix(
            src,
            output_file,
            str(input_file),
        )


def find_gse_folders(root: Path):
    """Return only folders beginning with GSE."""

    if root.is_dir() and root.name.upper().startswith("GSE"):
        return [root]

    return sorted(
        folder
        for folder in root.rglob("*")
        if folder.is_dir()
        and folder.name.upper().startswith("GSE")
    )


def extract_geo_matrices(
    input_path,
    output_dir=DEFAULT_OUTPUT_DIR,
):
    """Reusable pipeline function.

    Example
    -------
    from data.extract_series_matrix import extract_geo_matrices

    extract_geo_matrices("/path/to/GEO/downloads")
    """

    root = Path(input_path)
    output_dir = Path(output_dir)

    if not root.exists():
        raise FileNotFoundError(
            f"Input path not found: {root}"
        )

    gse_folders = find_gse_folders(root)

    if not gse_folders:
        raise RuntimeError(
            f"No GSE* folders found under {root}"
        )

    print(f"Found {len(gse_folders)} GSE folder(s).")

    created = []
    failures = []

    for gse_folder in gse_folders:

        print()
        print("=" * 60)
        print(f"Processing: {gse_folder.name}")
        print("=" * 60)

        matrix_files = sorted(
            file
            for file in gse_folder.rglob("*")
            if file.is_file()
            and "series_matrix" in file.name.lower()
            and (
                file.name.lower().endswith(".txt")
                or file.name.lower().endswith(".txt.gz")
            )
        )

        if not matrix_files:
            print("No series matrix found.")
            continue

        for matrix_file in matrix_files:

            try:
                output = process_matrix_file(
                    matrix_file,
                    output_dir,
                )

                created.append(output)

            except Exception as exc:

                failures.append(
                    (matrix_file, str(exc))
                )

                print()
                print(f"FAILED: {matrix_file}")
                print(f"Reason: {exc}")

    print()
    print("=" * 60)
    print("GEO extraction finished")
    print(f"GSE folders scanned: {len(gse_folders)}")
    print(f"Matrices created:    {len(created)}")
    print(f"Failures:            {len(failures)}")
    print(f"Output directory:    {output_dir}")
    print("=" * 60)

    return created


def main():

    parser = argparse.ArgumentParser(
        description=(
            "Extract GEO series matrices from GSE dataset folders."
        )
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        help="Parent directory containing GSE* folders",
    )

    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Output directory for intensity CSV files",
    )

    args = parser.parse_args()

    extract_geo_matrices(
        args.input_dir,
        args.output_dir,
    )


if __name__ == "__main__":
    main()