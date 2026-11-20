#!/usr/bin/env python3
"""Plot a stress-strain CSV produced by foam::write_csv.

Not part of the C++ build -- just a convenience script, since matplotlib
is a much faster way to eyeball a curve than writing a C++ plotting path.

Usage:
    python3 scripts/plot_stress_strain.py path/to/curve.csv
"""
import sys
import csv

import matplotlib.pyplot as plt


def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <csv_path>")
        sys.exit(1)

    path = sys.argv[1]
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        rows = [(float(x), float(y)) for x, y in reader]

    xs, ys = zip(*rows)
    plt.plot(xs, ys, marker="o", markersize=3)
    plt.xlabel(header[0])
    plt.ylabel(header[1])
    plt.title(path)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
