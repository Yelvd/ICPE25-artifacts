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

results_dir = "./extrapolation-results/"


w = walk(results_dir)
columns = ["jobid", "platform", "tasks", "benchmark", "tool"]
experiments = []
for (dirpath, dirnames, filenames) in w:
    if "meta.yaml" not in filenames:
        continue

    with open(f"{dirpath}/meta.yaml") as f:
        meta = yaml.load(f, Loader=yaml.FullLoader)
    
    jobid = meta["General"]["Id"]
    platform = meta["General"]["Platform"]
    tasks = meta["General"]["Tasks"]

    if meta["General"].get("Benchmark"):
        experiment = meta["General"]["Benchmark"]
    else:
        experiment = "fixed"
     
    if meta["General"].get("Monitor-tool"):
        tool = meta["General"]["Monitor-tool"]
    else:
        tool = meta["Instrumentation"]["Collection tool"]

    experiments.append([jobid, platform, tasks, experiment, tool])


df = pd.DataFrame(experiments, columns=columns)
# df = df.loc[((df['platform'] == "snellius-genoa") | (df['platform'] == "snellius-rome"))]
# df = df.loc[df['tool'] == "likwid"]
# print(df)

# for i in df['jobid']:
    # ear = False
    # d = f"./results/snellius/{i}"
    # files = [f for f in listdir(d) if isfile(join(d, f))]
    # for f in files:
        # if re.match("ear.csv", f):
            # ear = True

    # if not ear:
        # pass


df = df.groupby(["platform", "tasks", "benchmark", "tool"]).count()

with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    print(df)
