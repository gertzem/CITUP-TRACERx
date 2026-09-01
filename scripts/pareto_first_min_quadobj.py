import re
import argparse
import sys

import pandas as pd

cruk = re.compile(r"CRUKP?\d+(?:Tumor\d)?")

parser = argparse.ArgumentParser()

parser.add_argument(
    "--tau",
    type=float,
    default=1.05,
    help="Tau value (default: 1.05)",
)

parser.add_argument(
    "--epsilon",
    type=float,
    default=0.0001,
    help="Epsilon value (default: 0.0001)",
)

parser.add_argument(
    "files",
    nargs="*",
    help="Input files",
)

args = parser.parse_args()

epsilon = args.epsilon
tau = args.tau


optima = []
for filename in args.files:
    m = cruk.search(filename)
    if not m:
        print("Cannot read", filename)
        sys.exit(1)
    person = m.group()
    scores = pd.read_csv(filename, sep="\t", names=["migration", "tree_score"])
    best_tree_score = scores.tree_score.min()
    upper_limit = tau * (best_tree_score + epsilon)
    good_scores = scores.query("tree_score <= @upper_limit")
    i = good_scores.migration.argmin()

    optima.append((person, good_scores.migration.iloc[i], good_scores.tree_score.iloc[i]))

for person, migration, tree_score in sorted(optima):
    print(person, migration, tree_score, sep="\t")
