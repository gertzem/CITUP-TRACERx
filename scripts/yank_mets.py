import pandas as pd

persons = pd.read_csv("samplesheet.tsv", sep="\t").person
print("person\tmetastases")
for person in persons:
    sites = pd.read_csv(f"data/tracerx/{person}/sites.tsv", sep="\t")
    mets = sites.query('site != "p.lung"').site.to_list()
    print(person, ",".join(mets), sep="\t")
