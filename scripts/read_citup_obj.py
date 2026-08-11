import sys
import re

cruk = re.compile(r"CRUK\d+")

for filename in sys.argv[1:]:
    m = cruk.search(filename)
    if not m:
        print("Cannot read", filename)
        sys.exit(0)
    person = m.group()

    score = "NA"
    migrations = "NA"
    reading_score = False
    reading_migrations = True
    with open(filename) as fh:
        for line in fh:
            if reading_score:
                score = line.rstrip()
                reading_score = False
            elif reading_migrations:
                migrations = line.rstrip()
                reading_migrations = False
            elif line.startswith("TREE_SCORE_NO_REWARD"):
                reading_score = True
            elif line.startswith("MIKE_MIGRATION_SCORE"):
                reading_migrations = True
            else:
                reading_score = False
                reading_migrations = False

    print(person, migrations, score, sep="\t")
