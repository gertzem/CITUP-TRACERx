import sys
import re

cruk = re.compile(r"CRUKP?\d+(?:_Tumour\d)?")
refined_tree = re.compile(r"refined.tree")

def is_dag(children):
    # 0 = unvisited, 1 = visiting, 2 = finished
    state = {}

    def visit(node):
        if state.get(node, 0) == 1:
            # Found a back edge -> directed cycle
            return False
        if state.get(node, 0) == 2:
            # Already fully explored
            return True

        state[node] = 1  # currently visiting

        for child in children.get(node, []):
            if not visit(child):
                return False

        state[node] = 2  # fully explored
        return True

    # Include nodes that appear only as tails
    nodes = set(children)
    for tails in children.values():
        nodes.update(tails)

    return all(visit(node) for node in nodes)


for filename in sys.argv[1:]:
    locations = dict()
    location_file = re.sub(refined_tree, 'location.labeling', filename)
    with open(location_file) as fh:
        for line in fh:
            node, location = line.split()
            locations[node] = location
    m = cruk.search(filename)
    if not m:
        print("Cannot read", filename)
        sys.exit(0)
    person = m.group()

    migrations = 0
    comigrations = 0
    children = dict()
    with open(filename) as fh:
        for line in fh:
            head, tail, mig = line.split()
            mig = int(mig)
            if mig >= 0:
                comigrations = max(comigrations, mig + 1)
                migrations += 1
                head = re.sub(node, '', locations[head])
                tail = re.sub(node, '', locations[tail])
                children.setdefault(head, set()).add(tail)
    print(person, migrations, comigrations, is_dag(children), sep="\t")
