import sys
import re

cruk = re.compile(r"CRUK\d+")

for filename in sys.argv[1:]:
    m = cruk.search(filename)
    if not m:
        print("Cannot read", filename)
        sys.exit(0)
    person = m.group()

    migrations = 0
    comigrations = 0
    with open(filename) as fh:
        for line in fh:
            _, _, mig = line.split()
            mig = int(mig)
            if mig >= 0:
                comigrations = max(comigrations, mig + 1)
                migrations += 1
    print(person, migrations, comigrations, sep="\t")
