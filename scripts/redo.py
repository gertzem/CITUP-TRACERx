import fileinput

for line in fileinput.input():
    fields = line.split()
    print(f"results/migration_scan/meso/{fields[1]}/{fields[0]}")
