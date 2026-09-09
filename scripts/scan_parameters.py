import pandas as pd
import sys

def build_tests(center):
    start = max(center - 5, 1)
    finish = center + 5
    return ','.join([str(i) for i in range(start, finish+1)])

samplesheet = pd.read_csv(sys.argv[1], sep="\t")
samplesheet['orig_tests'] = samplesheet['tests']
samplesheet = samplesheet.drop(columns=['tests'])
penalty_optima = pd.read_csv(sys.argv[2], sep="\t",
                           names=['person', 'migrations', 'commigrations', 'score', 'tree_hash'])

penalty_optima['tests'] = penalty_optima['migrations'].apply(build_tests)
new_samplesheet = samplesheet.merge(penalty_optima[['person', 'tests']], how='left', on='person')
new_samplesheet.to_csv(sys.stdout, sep="\t", index=False)