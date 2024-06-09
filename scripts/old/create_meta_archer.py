import sys
import os
import yaml
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from bs4 import BeautifulSoup


results_dir = "./results-tmp/"
script = "scripts/old/generate-meta-data-yaml.py"

def which_experiment(experiment_dir):
    filename = glob.glob(f'{experiment_dir}/config*.xml')[0]
    # Reading the data inside the xml file to a variable under the name  data
    with open(filename, 'r') as file:
        data = file.read()

    # Passing the stored data inside the beautifulsoup parser
    bs_data = BeautifulSoup(data, 'xml')

    # Using find() to extract attributes of the first instance of the tag
    b_name = bs_data.find('nz')
    for test in b_name.children:
        nz = int(test)

    b_name = bs_data.find('ny')
    for test in b_name.children:
        ny = int(test)

    b_name = bs_data.find('tmax')
    for test in b_name.children:
        tmax = int(test)

    if ny == 200 and nz == 100:
        exp = "fixed"

    if ny == 200 and nz == 200:
        exp = "fixed-large"

    if ny == 100 and nz == 100:
        exp = "relative"

    if tmax == 0:
        return f"{exp}-empty"
    else:
        return f"{exp}"




columns = ["jobid", "platform", "tasks", "benchmark", "tool"]
experiments = []
for d in os.listdir(results_dir):
    exp_dir = f"{results_dir}/{d}/"
    files = [f for f in os.listdir(f"{exp_dir}/") if os.path.isfile(f"{exp_dir}/{f}")]
    
    meta = False
    perf = False
    for f in files:
        if "meta" in f:
            meta = True
        if "perf" in f:
            perf = True

    if meta:
        continue

    exp = which_experiment(exp_dir)
    jobid = d
    tool = "perf" if perf else "internal"

    cmd = f"python3 {script} -m {tool} -p archer_rome -i {jobid} -t 64 -n 1 -b {exp} -f {exp_dir}/meta.yaml ~/HemoCell-dev  ~/HemoCell-dev/Hemocell-Performance-Benchmarks/cube-benchmark"

    print(cmd)



# for (dirpath, dirnames, filenames) in w:
    # if "meta.yaml" not in filenames:
        # continue

    # with open(f"{dirpath}/meta.yaml") as f:
        # meta = yaml.load(f, Loader=yaml.FullLoader)
    
    # jobid = meta["General"]["Id"]
    # platform = meta["General"]["Platform"]
    # tasks = meta["General"]["Tasks"]

    # if meta["General"].get("Benchmark"):
        # experiment = meta["General"]["Benchmark"]
    # else:
        # experiment = "fixed"
     
    # if meta["General"].get("Monitor-tool"):
        # tool = meta["General"]["Monitor-tool"]
    # else:
        # tool = meta["Instrumentation"]["Collection tool"]

    # experiments.append([jobid, platform, tasks, experiment, tool])


# df = pd.DataFrame(experiments, columns=columns)

# df = df.groupby(["platform", "tasks", "benchmark", "tool"]).count()

# with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    # print(df)
