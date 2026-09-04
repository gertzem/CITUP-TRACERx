from pathlib import Path
import pandas as pd
import sys
import subprocess

samplesheet = pd.read_csv(sys.argv[1], sep="\t")
samplesheet = samplesheet[~samplesheet.argmin.isna()]

optima = sys.argv[2]

cohort_dir = Path(optima).parent.parent
cohort = cohort_dir.name
outdir = cohort_dir / "plots"

for _, r in samplesheet.iterrows():
    indir = cohort_dir / f"mig_{r['argmin']:03d}" / r.person

    subprocess.run(
        [
            "python",
            "../CITUP/plot_with_gene_names.py",
            indir,
            "--output-dir",
            outdir,
            "--sample-sites",
            f"data/{cohort}/{r.person}/sites.tsv",
            "--with-site-transitions",
            "--with_node_info",
            "--all-mike-dag-colorings",
        ],
        check=True,
    )

    subprocess.run(
        [
            "python",
            "../CITUP/plot_commigration_graphs.py",
            indir,
            "--output-dir",
            outdir,
            "--sites-file",
            f"data/{cohort}/{r.person}/sites.tsv",
        ],
        check=True,
    )
