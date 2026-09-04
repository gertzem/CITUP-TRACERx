import pandas as pd
import sys

samplesheet = pd.read_csv(sys.argv[1], sep='\t')
optima = pd.read_csv(sys.argv[2], sep="\t", names=["person", 'optimum', 'score'])
optima.merge(samplesheet, on="person").query('optimum != argmin')[['person', 'optimum', 'argmin']].to_csv(sys.stdout,sep="\t", index=False)
