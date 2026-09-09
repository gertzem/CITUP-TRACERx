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
#     display_name: Python (base)
#     language: python
#     name: base
# ---

import pandas as pd

us = pd.read_csv(
    "results/migration_scan/tracerx/pareto/optima.tsv",
    names=["person", "migrations", "comigrations", "quad", "tree"],
    sep="\t",
)

them = pd.read_csv(
    "results/fixed_trees/optima.tsv",
    names=["person", "migrations", "comigrations", "quad", "tree"],
    sep="\t",
)

both = us.merge(them, on="person", suffixes=("_us", "_them"))

mach2 = pd.read_csv(
    "results/summary/mach2_migrations.tsv",
    sep="\t",
    names=["person", "migrations_mach2", "comigrations_mach2", "is_dag"],
)

cmp = both.merge(mach2, on="person").set_index("person", drop=False)


def migrations_cmp(a, b):
    if a < b:
        return "migrations_better"
    elif a > b:
        return "migrations_worse"
    else:
        return "migrations_same"


def mach2_comigrations_cmp(a, b):
    if a < b:
        return "0mach2_worse_cm"
    elif a > b:
        return "2mach2_better_cm"
    else:
        return "1mach2_same_cm"


def tree_score_cmp(a, b):
    if a < b:
        return "score_better"
    elif a > b:
        return "score_worse"
    else:
        return "score_same"


def tree_score_relcmp(a, b):
    if a < 0.95 * b:
        return "score_better"
    elif a > 1.05 * b:
        return "score_worse"
    else:
        return "score_same"


def mach2_cmp(a, b):
    if a < b:
        return "0mach2_worse"
    elif a > b:
        return "2mach2_better"
    else:
        return "1mach2_same"


cmp["migrations_cmp"] = cmp.apply(
    lambda r: migrations_cmp(r["migrations_us"], r["migrations_them"]), axis=1
)

cmp["fixed_mig_cmp"] = cmp.apply(
    lambda r: migrations_cmp(r["migrations_them"], r["migrations_mach2"]),
    axis=1,
)

cmp["mach2_cmp"] = cmp.apply(
    lambda r: mach2_cmp(r["migrations_us"], r["migrations_mach2"]), axis=1
)

cmp["mach2_comigrations_cmp"] = cmp.apply(
    lambda r: mach2_comigrations_cmp(
        r["comigrations_us"], r["comigrations_mach2"]
    ),
    axis=1,
)

cmp["tree_score_cmp"] = cmp.apply(
    lambda r: tree_score_cmp(r["quad_us"], r["quad_them"]), axis=1
)

cmp["tree_score_relcmp"] = cmp.apply(
    lambda r: tree_score_relcmp(r["quad_us"], r["quad_them"]), axis=1
)

# ## Absolute score comparisons

print(pd.crosstab(cmp["migrations_cmp"], cmp["tree_score_cmp"]).to_markdown())

print(pd.crosstab(cmp["mach2_cmp"], cmp["tree_score_cmp"]).to_markdown())

# ## Relative Score comparisions

print(
    pd.crosstab(cmp["migrations_cmp"], cmp["tree_score_relcmp"]).to_markdown()
)

print(pd.crosstab(cmp["mach2_cmp"], cmp["tree_score_relcmp"]).to_markdown())

# ## There are no cases on the fixed CONIPHER tree that mach2 gets fewer migrations

cmp.query('fixed_mig_cmp == "mach2_better"')

print(
    cmp.query('tree_score_relcmp == "score_worse"')[
        [
            "person",
            "migrations_us",
            "quad_us",
            "migrations_them",
            "quad_them",
            "migrations_mach2",
        ]
    ].to_markdown(index=False)
)

cmp.query("tree_us == tree_them").shape[0]

cmp.query("migrations_us == migrations_them & quad_us == quad_them").shape[0]

cmp["same_tree"] = cmp["tree_us"] == cmp["tree_them"]

print(
    cmp.query('mach2_cmp == "0mach2_worse" & tree_score_cmp == "score_same"')[
        ["person", "migrations_us", "migrations_mach2", "same_tree"]
    ].to_markdown(index=False)
)

print(
    cmp.query(
        'mach2_comigrations_cmp == "2mach2_better_cm" & mach2_cmp == "1mach2_same" &tree_score_cmp == "score_same"'
    )[
        [
            "person",
            "migrations_us",
            "comigrations_us",
            "comigrations_mach2",
            "same_tree",
        ]
    ].to_markdown(index=False)
)

mach2.is_dag.value_counts()
