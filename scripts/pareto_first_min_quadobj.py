import re
import sys

import pandas as pd

cruk = re.compile(r"CRUK\d+")
epsilon = 1e-4
tau = 1.05

optima = []
for filename in sys.argv[1:]:
    m = cruk.search(filename)
    if not m:
        print("Cannot read", filename)
        sys.exit(0)
    person = m.group()
    scores = pd.read_csv(filename, sep="\t", names=["migration", "tree_score"])
    best_tree_score = scores.tree_score.min()
    upper_limit = tau * (best_tree_score + epsilon)
    good_scores = scores.query("tree_score <= @upper_limit")
    i = good_scores.migration.argmin()

    optima.append((person, good_scores.migration.iloc[i], good_scores.tree_score.iloc[i]))

for person, migration, tree_score in sorted(optima):
    print(person, migration, tree_score, sep="\t")
