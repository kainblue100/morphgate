"""Command line for a HER2ST leave-one-patient run."""

from __future__ import annotations

import argparse
from pathlib import Path

from morphgate.fit import fit_her2st


def main() -> None:
    parser = argparse.ArgumentParser(description="Predict morphology-linked spatial expression from H&E.")
    parser.add_argument("--root", type=Path, required=True, help="Directory that contains the unpacked HER2ST files.")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    fit_her2st(args.root, args.out)


if __name__ == "__main__":
    main()
