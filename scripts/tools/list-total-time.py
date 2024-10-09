import sys
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from os import walk
from os import listdir
from os.path import isfile, join
from os import walk
import re

results_dir = "./results/"


w = walk(results_dir)
columns = ["jobid", "platform", "tasks", "benchmark", "tool"]
total = 0
for (dirpath, dirnames, filenames) in w:

    filename = None
    for file in filenames:
        if ".out" in file:
            filename = file
            break

    if filename:
        with open(f"{dirpath}/{filename}") as f:
            for line in f:
                if "HemoCell: " in line:
                    total = total + float(line.split()[-1])
                    break

print(f"Total seconds: {total}")
print(f"Total min: {total / 60}")
print(f"Total hour: {total / 60 / 60}")
