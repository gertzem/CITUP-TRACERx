import argparse
import csv
import gzip
import pickle
from collections import Counter
from typing import NamedTuple
import numpy as np
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--type", required=True, choices=["citupm", "mach2", "metient"])
parser.add_argument("ground_truth")
parser.add_argument("inferred", nargs="*")

args = parser.parse_args()


class Edge(NamedTuple):
    source: str
    target: str
    copy: int


def read_ground_truth(path):
    counts = Counter()
    edges = set()
    with open(path, newline="") as file:
        for line in file:
            source, target = line.split()
            counts[source, target] += 1
            edges.add(Edge(source, target, counts[source, target]))
    return edges


def read_mach2(path: str) -> set[Edge]:
    edges = set()
    with open(path, newline="") as file:
        for source, target, multiplicity in csv.reader(file, delimiter="\t"):
            for copy in range(1, int(multiplicity) + 1):
                edges.add(Edge(source, target, copy))
    return edges


def read_citupm(path: str) -> set[Edge]:
    counts = Counter()
    edges = set()
    with open(path, newline="") as file:
        for row in csv.DictReader(file, delimiter="\t"):
            if row["type"] != "migration":
                continue
            source = row["from_site"]
            target = row["to_site"]
            counts[source, target] += 1
            edges.add(Edge(source, target, counts[source, target]))
    return edges


def read_metient(path: str) -> set[Edge]:
    import metient as met
    with gzip.open(path, "rb") as fh:
        pkl = pickle.load(fh)
    sites = pkl['anatomical_sites']
    V = pkl["node_labels"][0]
    A = met.adjacency_matrix_from_parents(pkl["parents"][0])
    G = met.migration_graph(V, A)
    edges = set()
    for i, source in enumerate(sites):
        for j, target in enumerate(sites):
            if G[i][j] != 0:
                edges.add(Edge(source, target, int(G[i][j])))
    return edges

case = Path(args.ground_truth).parent.name
ground_truth = read_ground_truth(args.ground_truth)
scores = []
for path in args.inferred:
    if args.type == 'citupm':
        inferred = read_citupm(path)
    elif args.type == 'mach2':
        inferred = read_mach2(path)
    else:
        inferred = read_metient(path)
    common = len(ground_truth.intersection(inferred))
    recall = common / len(ground_truth)
    precision = common / len(inferred)
    if precision + recall == 0:
        f1 = 0
    else:
        f1 = 2 * (precision * recall) / (precision + recall)
    scores.append(f1)
    # print(precision, recall, f1, sep="\t")
print(case, np.median(scores), np.median(precision), np.median(recall), sep="\t")