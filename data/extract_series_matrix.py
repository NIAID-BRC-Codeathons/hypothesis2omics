"""Extract the expression matrix from a GEO series-matrix TXT file.

The script extracts everything between:

    !series_matrix_table_begin
    !series_matrix_table_end

The marker lines themselves are excluded.

Output is saved as CSV in:
    data/galaxy_file_input/

The GEO first-column name ID_REF is renamed to probe_id.
"""

import argparse
import csv
from pathlib import Path


DEFAULT_OUTPUT_DIR = Path("data/galaxy_file_input")


def extract_matrix(input_file: Path, output_file: Path) -> None:
    inside_matrix = False
    found_end = False
    rows_written = 0

    output_file.parent.mkdir(parents=True, exist_ok=True)

    with input_file.open(
        "r",
        encoding="utf-8",
        errors="replace",
        newline="",
    ) as src, output_file.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as dst:

        writer = csv.writer(dst)

        for line in src:
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

            if rows_written == 0 and row and row[0] == "ID_REF":
                row[0] = "probe_id"

            writer.writerow(row)
            rows_written += 1

    if not inside_matrix:
        raise RuntimeError(
            "Could not find !series_matrix_table_begin in the input file."
        )

    if not found_end:
        raise RuntimeError(
            "Could not find !series_matrix_table_end in the input file."
        )

    if rows_written == 0:
        raise RuntimeError("No matrix rows were extracted.")

    print("Matrix extraction complete")
    print(f"Input:  {input_file}")
    print(f"Output: {output_file}")
    print(f"Rows written including header: {rows_written}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract a GEO series matrix for Galaxy input."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Path to GEO series_matrix.txt file",
    )

    parser.add_argument(
        "--output-name",
        required=True,
        help="Output CSV filename",
    )

    args = parser.parse_args()

    input_file = Path(args.input)

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    output_file = DEFAULT_OUTPUT_DIR / args.output_name

    extract_matrix(input_file, output_file)


if __name__ == "__main__":
    main()