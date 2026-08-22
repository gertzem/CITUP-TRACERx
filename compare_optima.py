# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.19.3
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

import pandas as pd

us = pd.read_csv('results/citup_scan_mets/conipher/pareto/optima.tsv',
                 names=['person', 'migrations', 'quad', 'tree'], sep="\t")

them = pd.read_csv('results/citup2_mets/conipher/optima.tsv',
                 names=['person', 'migrations', 'quad', 'tree'], sep="\t")

both = us.merge(them, on='person', suffixes=("_us", "_them"))

mach2 = pd.read_csv('results/summary/mach2_migrations.tsv', sep="\t", names=["person", "migrations_mach2", "comigrations_mach2"])

mach2

cmp = both.merge(mach2, on='person').set_index('person', drop=False)


def migrations_cmp(a, b):
    if a < b:
        return 'migrations_better'
    elif a > b:
        return 'migrations_worse'
    else:
        return 'migrations_same'


def mach2_cmp(a, b):
    if a < b:
        return '0mach2_worse'
    elif a > b:
        return '2mach2_better'
    else:
        return '1mach2_same'


def tree_score_cmp(a, b):
    if a < b:
        return 'score_better'
    elif a > b:
        return 'score_worse'
    else:
        return 'score_same'


def tree_score_relcmp(a, b):
    if a < 0.95 * b:
        return 'score_better'
    elif a > 1.05 * b:
        return 'score_worse'
    else:
        return 'score_same'


cmp['migrations_cmp'] = cmp.apply(lambda r: migrations_cmp(r['migrations_us'], r['migrations_them']), axis=1)

cmp['mach2_cmp'] = cmp.apply(lambda r: mach2_cmp(r['migrations_us'], r['migrations_mach2']), axis=1)

cmp['tree_score_cmp'] = cmp.apply(lambda r: tree_score_cmp(r['quad_us'], r['quad_them']), axis=1)

cmp['tree_score_relcmp'] = cmp.apply(lambda r: tree_score_relcmp(r['quad_us'], r['quad_them']), axis=1)

# Absolute score comparisons

pd.crosstab(cmp['migrations_cmp'], cmp['tree_score_cmp'])

pd.crosstab(cmp['mach2_cmp'], cmp['tree_score_cmp'])

pd.crosstab(cmp['migrations_cmp'], cmp['tree_score_relcmp'])

pd.crosstab(cmp['mach2_cmp'], cmp['tree_score_relcmp'])

# Relative score commparisons

cmp.query('tree_score_relcmp == "score_worse"')

cmp.query('migrations_cmp == "migrations_worse" & mach2_cmp == "2mach2_better"')

cmp.query('migrations_cmp == "migrations_better"')

cmp.query('tree_us == tree_them').shape[0]

cmp.query('tree_us == tree_them & quad_us == quad_them').shape[0]


