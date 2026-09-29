import pandas as pd
import sys

samplesheet = pd.read_csv('samplesheet-tracerx.tsv', sep="\t")
for person in samplesheet.person:

    pareto = pd.read_csv(f"results/migration_scan/tracerx/pareto/{person}/pareto.txt",
                         sep="\t",
                         names=["person", "migrations", "comigrations", "tree_score", "tree_hash"])
    start = pareto.migrations.min()
    finish = pareto.migrations.max()
    migrations = set(pareto.migrations)
    for i in range(start, finish+1):
        if i not in migrations:
            print(person, i, sep="\t")