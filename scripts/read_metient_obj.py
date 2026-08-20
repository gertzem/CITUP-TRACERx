import sys
import re
import gzip
import pickle

cruk = re.compile(r"CRUK\d+")

for filename in sys.argv[1:]:
    m = cruk.search(filename)
    if not m:
        print("Cannot read", filename)
        sys.exit(0)
    person = m.group()

    with gzip.open(filename, 'rb') as fh:
        met = pickle.load(fh)
        migrations = met['loss_info'][0][
            'migration_number'
        ]
        comigrations = met['loss_info'][0][
            'comigration_number'
        ]
    print(person, migrations, comigrations, sep="\t")
