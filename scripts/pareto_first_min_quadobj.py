#!/usr/bin/env python3

import sys
from pathlib import Path

for filename in sys.argv[1:]:
    min_value = float("inf")
    argmin = None

    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            i, m = line.split()
            m = float(m)  # Use int(m) if the second column is always an integer

            if m < min_value:
                min_value = m
                argmin = i

    print(Path(filename).parent, argmin, min_value, sep="\t")
