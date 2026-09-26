import pandas as pd
import sys
import sys

cohort = sys.argv[1]

samplesheet = pd.read_csv(f'samplesheet-{cohort}.tsv', sep="\t")
sequential = pd.read_csv(f'results/migration_trees/{cohort}/optima.tsv', sep="\t",
                        names=['person', 'migrations', 'comigrations', 'tree_score', 'tree_hash'])
sequential = sequential.set_index('person', drop=False)
print('tests')
for person in samplesheet.person:
    sites = pd.read_csv(f'data/{cohort}/{person}/sites.tsv', sep="\t")
    nsites = sites.site.nunique()
    print(','.join([str(i) for i in range(nsites - 1, sequential.loc[person, 'migrations']+1)]))