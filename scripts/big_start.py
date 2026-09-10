import pandas as pd
import sys

samplesheet = pd.read_csv(sys.argv[1], sep="\t").set_index('person', drop=False)
smallest_test = samplesheet.tests.apply(lambda x: min(map(int, x.split(","))))
for person, smallest in smallest_test.items():
    sites = pd.read_csv(f"data/tracerx/{person}/sites.tsv",sep="\t")
    actual_min = sites.site.nunique() - 1
    if smallest > actual_min:
        for i in range(actual_min, smallest):
            print(f"results/migration_scan/tracerx/mig_{i:03d}/{person}")