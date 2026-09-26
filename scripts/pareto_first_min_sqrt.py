import re
import argparse
import sys
import math
import numpy as np

import pandas as pd

cruk = re.compile(r"P\d+a?")

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
    "baseline_optima"
)

parser.add_argument(
    "files",
    nargs="*",
    help="Input files",
)

args = parser.parse_args()

epsilon = args.epsilon
tau = args.tau

baseline_optima = pd.read_csv(args.baseline_optima, sep="\t", names=["person", "migrations", "commigrations", "tree_score", "tree_hash"])
baseline_optima = baseline_optima.set_index('person', drop=False)

optima = []
for filename in args.files:
    m = cruk.search(filename)
    if not m:
        print("Cannot read", filename)
        sys.exit(1)
    person = m.group()
    baseline = baseline_optima.loc[person]
    scores = pd.read_csv(filename, sep="\t", names=["person", "migrations", "commigrations", "tree_score", "tree_hash"])
    scores = scores.set_index("migrations", drop=False)
    scores = scores[scores.migrations <= baseline.migrations]
    if baseline.migrations in scores.index:
        if baseline.tree_score < scores.loc[baseline.migrations, 'tree_score']:
            scores.loc[baseline.migrations, :] = baseline
    else:
        scores.loc[baseline.migrations, :] = baseline

    best_tree_score = scores.tree_score.min()
    upper_limit = tau * (math.sqrt(best_tree_score) + epsilon)
    good_scores = scores[np.sqrt(scores['tree_score']) <= upper_limit]
    i = good_scores.migrations.argmin()

    optima.append(good_scores.iloc[i])

optima = pd.DataFrame(optima).sort_values('person')
optima.to_csv(sys.stdout, sep="\t", index=False, header=False)
