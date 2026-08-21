import pandas as pd
import json
import subprocess

samples = pd.DataFrame(json.loads(snakemake.params.samples))

arglist = []
for _, sample in samples.iterrows():
    arglist.append(f"results/citup_scan_mets/conipher/mig_{sample['argmin']:03}/{sample.person}/{sample.person}.mig_{sample['argmin']:03}.tree_0.txt")

with open(snakemake.output[0], 'w') as fh:
    subprocess.run(
        ['python', 'scripts/read_citup_obj.py', *arglist],
        stdout=fh,
        check=True
    )