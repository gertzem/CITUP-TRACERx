import re
import argparse
import sys
import math
import numpy as np

import pandas as pd

cruk = re.compile(r"CRUKP?\d+(?:_Tumour\d)?")

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
    scores = pd.read_csv(filename, sep="\t", names=["person", "migration", "commigrations", "tree_score", "tree_hash"])
    best_tree_score = scores.tree_score.min()
    upper_limit = tau * (math.sqrt(best_tree_score) + epsilon)
    good_scores = scores[np.sqrt(scores['tree_score']) <= upper_limit]
    i = good_scores.migration.argmin()

    optima.append(good_scores.iloc[i])

optima = pd.DataFrame(optima).sort_values('person')
optima.to_csv(sys.stdout, sep="\t", index=False, header=False)
