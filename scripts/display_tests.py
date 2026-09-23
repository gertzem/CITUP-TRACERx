import pandas as pd
import sys

samplesheet = pd.read_csv('samplesheet-tracerx.tsv', sep="\t")
sequential = pd.read_csv('results/migration_trees/tracerx/optima.tsv', sep="\t",
                        names=['person', 'migrations', 'comigrations', 'tree_score', 'tree_hash'])
sequential = sequential.set_index('person', drop=False)
print('display_tests')
for person in samplesheet.person:
    sites = pd.read_csv(f'data/tracerx/{person}/sites.tsv', sep="\t")
    print(','.join([str(i) for i in range(sequential.loc[person, 'migrations']+1, sequential.loc[person, 'migrations']+6)]))