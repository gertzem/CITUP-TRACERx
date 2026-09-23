import pandas as pd
from pathlib import Path

sequential = pd.read_csv('results/migration_trees/tracerx/optima.tsv', sep="\t",
                        names=['person', 'migrations', 'comigrations', 'tree_score', 'tree_hash'])

pareto = []

for filename in Path('results/migration_scan/tracerx/pareto').rglob('pareto.txt'):
    pareto.append(pd.read_csv(filename, sep="\t",
                        names=['person', 'migrations', 'comigrations', 'tree_score', 'tree_hash']))

pareto = pd.concat(pareto)

commute = sequential[['person', 'migrations', 'tree_score']].merge(pareto[['person', 'migrations', 'tree_score']], how='left', on=['person', 'migrations'])

print('Non-commuting:')
print(commute.query('tree_score_y.isna() | (tree_score_x != tree_score_y)').to_markdown())
