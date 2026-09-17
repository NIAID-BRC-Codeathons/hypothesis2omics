"""Extract GEO series-matrix expression tables for Galaxy input.

Supports:
1. A single GEO series-matrix TXT file with --input
2. All *_series_matrix.txt files in a folder with --input-dir

Outputs are written to:
    data/galaxy_file_input/

Example output:
    GSE13699_GPL6104_intensities.csv
"""

import argparse
import csv
from pathlib import Path


DEFAULT_OUTPUT_DIR = Path("data/galaxy_file_input")


def make_output_name(input_file: Path) -> str:
    """Create a Galaxy-friendly output filename."""

    name = input_file.name

    if name.endswith("_series_matrix.txt"):
        name = name[: -len("_series_matrix.txt")]

    name = name.replace("-", "_")

    return f"{name}_intensities.csv"


def extract_matrix(input_file: Path, output_file: Path) -> None:
    """Extract the table between GEO series-matrix begin/end markers."""

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
            f"Could not find !series_matrix_table_begin in {input_file}"
        )

    if not found_end:
        raise RuntimeError(
            f"Could not find !series_matrix_table_end in {input_file}"
        )

    if rows_written == 0:
        raise RuntimeError(
            f"No matrix rows were extracted from {input_file}"
        )

    print()
    print("Matrix extraction complete")
    print(f"Input:  {input_file}")
    print(f"Output: {output_file}")
    print(f"Rows written including header: {rows_written}")


def process_single_file(
    input_file: Path,
    output_dir: Path,
    output_name: str | None = None,
) -> None:

    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")

    if output_name is None:
        output_name = make_output_name(input_file)

    output_file = output_dir / output_name

    extract_matrix(input_file, output_file)


def process_directory(input_dir: Path, output_dir: Path) -> None:
    """Process every GEO series-matrix TXT file in a directory."""

    if not input_dir.exists():
        raise FileNotFoundError(
            f"Input directory not found: {input_dir}"
        )

    matrix_files = sorted(
        input_dir.glob("*_series_matrix.txt")
    )

    if not matrix_files:
        raise RuntimeError(
            f"No *_series_matrix.txt files found in {input_dir}"
        )

    print(f"Found {len(matrix_files)} series matrix file(s).")

    successful = 0
    failed = 0

    for input_file in matrix_files:
        try:
            process_single_file(
                input_file=input_file,
                output_dir=output_dir,
            )
            successful += 1

        except Exception as exc:
            failed += 1
            print()
            print(f"FAILED: {input_file}")
            print(f"Reason: {exc}")

    print()
    print("Batch extraction finished")
    print(f"Successful: {successful}")
    print(f"Failed:     {failed}")
    print(f"Output directory: {output_dir}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Extract GEO series matrices for Galaxy input."
    )

    group = parser.add_mutually_exclusive_group(required=True)

    group.add_argument(
        "--input",
        help="Single GEO series_matrix.txt file",
    )

    group.add_argument(
        "--input-dir",
        help="Directory containing *_series_matrix.txt files",
    )

    parser.add_argument(
        "--output-name",
        help="Optional output filename when processing one file",
    )

    parser.add_argument(
        "--output-dir",
        default=str(DEFAULT_OUTPUT_DIR),
        help="Output directory",
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)

    if args.input:
        process_single_file(
            input_file=Path(args.input),
            output_dir=output_dir,
            output_name=args.output_name,
        )

    else:
        process_directory(
            input_dir=Path(args.input_dir),
            output_dir=output_dir,
        )


if __name__ == "__main__":
    main()