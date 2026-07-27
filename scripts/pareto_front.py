import sys

scores = []
for filename in sys.argv[1:]:
    with open(filename) as fh:
        for i, line in enumerate(fh):
            if i == 1:
                fields = line.split()
                scores.append((int(fields[4]), fields[2]))

scores.sort()
for mig, score in scores:
    print(mig, score, sep="\t")
