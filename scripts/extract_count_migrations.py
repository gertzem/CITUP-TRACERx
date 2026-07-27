import pandas as pd
import sys
from pathlib import Path

for file in sys.argv[1:]:
    basedir = Path(file)
    base = basedir.name

    with open(basedir / "trees" / f"{base}.tree_0.txt") as fh:
        lines = list(fh)
    imigrations = lines.index("MIKE_MIGRATION_SCORE\n")
    mid = int(lines[imigrations + 1])
    start = max(mid - 3, 1)
    finish = mid + 3 + 1
    tests = ",".join([str(i) for i in range(start, finish)])
    results.append((base, mid, tests))

results = pd.DataFrame(results, columns=["person", "migrations", "tests"])
samples = samples.merge(results, on="person")
samples.to_csv(sys.stdout, sep="\t", index=False)
