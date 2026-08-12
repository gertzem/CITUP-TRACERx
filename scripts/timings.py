import sys

for filename in sys.argv[1:]:
    with open(filename) as fh:
        for line in fh:
            if line.startswith('Explored'):
                fields = line.split()
                components = filename.split('/')
                migrations = components[5].split('.')[0]
                print(components[3] + ":" + migrations[3:], fields[7], sep='\t')
