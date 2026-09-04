from pathlib import Path
import pandas as pd
import sys
import subprocess



for indir in sys.argv[1:]:
    path = indir.split('/')
    if len(path) < 4:
        print(indir, 'is not a valid plotting path')
    cohort = path[2]
    person = path[3]
    outdir = Path(indir).parent / 'plots'

    subprocess.run(
        [
            "python",
            "../CITUP/plot_with_gene_names.py",
            indir,
            "--output-dir",
            outdir,
            "--sample-sites",
            f"data/{cohort}/{person}/sites.tsv",
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
            f"data/{cohort}/{person}/sites.tsv",
        ],
        check=True,
    )
